import numpy as np, scipy as sp, sys, os, myvariant, mygene, json, multiprocessing, pysam, subprocess
import re, textwrap, networkx as nx, community, requests, argparse, uuid, pandas as pd, glob, shutil
import pickle
from joblib import Parallel, delayed
from networkx.utils import open_file, make_str
from Bio import SeqIO, Seq, SeqRecord
from Bio.SeqRecord import SeqRecord
from collections import Counter
from networkx.readwrite import json_graph
from pyliftover import LiftOver
import plotly.graph_objects as go
from collections import OrderedDict
import plotly.io as pio

#Function structure:
#Helper / resource functions
#SNP/CNV functions
#Open regions functions
#Insulated Neighborhood functions
#Whole genome functions
#Graphing functions

#Iterate over input DNA-interaction graph and fetch out all nodes which are "targets" depending on input parameters
def get_target_nodes(G,args,cell_type):
    target_nodes = list()
    if args["top_gene_fpkm"] is not None:
        target_nodes = [node[0] for node in G.nodes(data=True) if "fpkm" in node[1] and node[1]["fpkm"] > float(args["top_gene_fpkm"])]
    elif args["marker_region"] is not None:
        marker_arr = re.split(r"[-:]",args["marker_region"])
        chr_marker = marker_arr[0]
        marker_start = int(marker_arr[1])
        marker_end = int(marker_arr[2])
        for node in list(G.nodes(data=True)):
            arr = re.split(r"[-:]",node[0])
            node_chr = arr[0]
            node_start = int(arr[1])
            node_end = int(arr[2])
            if "fpkm" in node[1] and node_chr == chr_marker and node_start >= marker_start and node_end <= marker_end:
                target_nodes.append(node[0])
    elif args["insulated_neigh"] is not None:
        target_names = args["insulated_neigh"].strip(",").split(",")
        for node in list(G.nodes(data=True)):
            if "name" in node[1]:
                node_names = set(node[1]["name"].split(","))
                if len(node_names.intersection(target_names)) > 0:
                    target_nodes.append(node[0])
    else:
        if cell_type == "lung":
            source_gene_names = "MYC,ASCL1,MYCN,MYCL,P53"
        elif cell_type == "lung-mesc":
            source_gene_names = "myc,ascl1,mycn,mycl,p53"
        else:
            raise NameError("Cant find that cell type, try: lung, lung-mESC")
        default_gene_names = [x.upper() for x in source_gene_names.split(",")]
        for node in list(G.nodes(data=True)):
            if "name" in node[1]:
                node_names = node[1]["name"].split(",")
                if len(set(node_names).intersection(default_gene_names)) > 0:
                    target_nodes.append(node[0])
    return target_nodes

#Function to take in network edge and spit out BED format edge
def edge2_strBed(edge):
    arr_l = re.split(r"[-:]",edge[0])
    arr_r = re.split(r"[-:]",edge[1])
    if int(arr_l[1]) > int(arr_r[1]):
        temp = arr_r
        arr_r = arr_l
        arr_l = temp
    len_l = abs(int(arr_l[1])-int(arr_l[2]))
    len_r = abs(int(arr_r[1])-int(arr_r[2]))
    if arr_l[0] == arr_r[0]:
        out_str = arr_l[0] + "\t" + arr_l[1] + "\t" + arr_r[2] + "\tInteraction\t1000\t+\t" + \
        arr_l[1] + "\t" + arr_l[1] + "\t0,0,0\t2\t" + str(len_l) + "," + str(len_r) + \
        "\t0," + str(abs(int(arr_l[1])-int(arr_r[1])))
    else:
        out_str = ""
    return out_str

#Function to take in network edge and spit out Washu format edge
def edge2_strWashu(edge):
    arr_l = re.split(r"[-:]",edge[0])
    arr_r = re.split(r"[-:]",edge[1])
    if int(arr_l[1]) > int(arr_r[1]):
        temp = arr_r
        arr_r = arr_l
        arr_l = temp
    mid_l = int(arr_l[1]) + (abs(int(arr_l[1])-int(arr_l[2]))/2.0)
    mid_r = int(arr_r[1]) + (abs(int(arr_r[1])-int(arr_r[2]))/2.0)
    if arr_l[0] == arr_r[0]:
        out_str = arr_l[0] + ":" + str(mid_l) + "-" + str(mid_l+1) + "\t" + arr_l[0] + \
        ":" + str(mid_r) + "-" + str(mid_r+1)  + "\t" + str(dna_int_graph.get_edge_data(*edge)["weight"])
    else:
        out_str = ""
    return out_str

#Make BED file of all nodes w/ in a given graph G
def uniq_nodes(G):
    #put all nodes in bed, sort, and uniq in all_nodes_sorted.bed
    with open("all_nodes.bed","w") as all_nodes_bed:
        for node in G.nodes():
            arr = re.split(r"[-:]",node)
            all_nodes_bed.write(arr[0] + "\t" + str(arr[1]) + "\t" + str(arr[2]) + "\n")
    all_nodes_bed.close()
    os.system("sort -u all_nodes.bed | sort-bed - > all_nodes_sorted.bed")
    os.remove("all_nodes.bed")
    return G

#Overlap bed file w/ graph and return list of overlapping nodes
def overlap_graph_bed(G, bed):
    os.system("sort-bed " + bed + " | bedops -e -1 " + "all_nodes_sorted.bed - > " + "overlapped_nodes.bed")
    overlapped_nodes = set()
    with open("overlapped_nodes.bed","r") as over_bed:
        for node in over_bed:
            arr = node.strip().split()
            overlapped_nodes.add(str(arr[0]) + ":" + str(arr[1]) + "-" + str(arr[2]))
    os.remove("overlapped_nodes.bed")
    return overlapped_nodes

#Take graph G and label all nodes that are close enough to TSSs as genes and add relevant FPKM
def label_gene_nodes(G, in_fpkm, tss_ensembl):
    #overlap TSSs file against expression file, match node to TSS, label node w/ data 
    tran2fpkm = dict()
    tran2hgnc = dict()
    tran2chr = dict()
    tran2start = dict()

    with open(tss_ensembl, "r") as tran_starts:
        #Remove header line
        tran_starts.readline()
        for tss in tran_starts:
            if not tss.startswith("#"):
                arr = tss.strip().split()
                tran2chr[arr[3]] = "chr"+arr[5]
                tran2start[arr[3]] = int(arr[4])
                tran2hgnc[arr[3]] = arr[0]

    with open(in_fpkm, "r") as tran_fpkm:
        for tss in tran_fpkm:
            arr = tss.strip().split()
            tran = arr[0].split(".")[0]
            tran2fpkm[tran] = arr[6]

    bed_out = "fpkm2tss2gene.bed"
    count_skipped = 0
    with open(bed_out, "w") as fpkm_out:
        for tran in tran2fpkm.keys():
            if (tran in tran2hgnc.keys()):
                out_line = tran2chr[tran] + "\t" + str(tran2start[tran]) + "\t" + str(tran2start[tran]+1) + "\t" + tran2hgnc[tran] + "\t" + tran2fpkm[tran] + "\n"
                fpkm_out.write(out_line)
            else:
                count_skipped += 1
    
    tss_slop = 50000
    os.system("sort-bed fpkm2tss2gene.bed | closest-features --dist --closest - all_nodes_sorted.bed > closest_nodes.bed")
    closest_nodes = list()
    closest_names = list()
    #Assign each TSS to the closest cohesin peak less than 50kb away
    with open("closest_nodes.bed") as closest_nodes:
        for node in closest_nodes:
            arr = node.strip().split("|")
            if (arr[1] != "NA") and (abs(int(arr[2])) <= tss_slop) :
                node_arr = arr[1].split()
                gene_tss_arr = arr[0].split()
                node_name = (str(node_arr[0]) + ":" + str(node_arr[1]) + "-" + str(node_arr[2]))
                if "name" not in G.nodes[node_name]:
                    G.nodes[node_name]["name"] = set()
                    G.nodes[node_name]["fpkm"] = 0.0
                G.nodes[node_name]["gene"] = True
                tss_list = set((gene_tss_arr[3].upper()).split(","))
                G.nodes[node_name]["name"].update(tss_list)
                #Fetch gene annotation from mygene.info
#                mg = mygene.MyGeneInfo()
#                target_gene_out = mg.query("symbol:"+gene_tss_arr[3].upper(), fields="ensembl,exac", species="human", fetch_all=True)
                
                G.nodes[node_name]["fpkm"] += float(gene_tss_arr[4])
    #join gene names
    for node in list(G.nodes(data=True)):
        if ("name" in node[1]):
            G.nodes[node[0]]["name"] = ",".join(node[1]["name"])
    os.remove("fpkm2tss2gene.bed")
    os.remove("closest_nodes.bed")
    return G

def annotate_IN(target_node_name, dendogram_com, level, dna_int_graph, add_open, add_snp, homo_gen,
                open_peaks_file, picked_tissue, ctcf2bounds, comm_open_out, comm_var_out, comm_nodes_out, comm_edges_out, comm_edges_out_bed, out_dir):    
    #fetch Insulated Neighborhood and annotate
    if level <= 1:
        comm_out_path = out_dir + run_name + "_" + dna_int_graph.nodes[target_node_name]["in_name"].split(":")[0][0:80] + "_lvl"+str(level+1)
    else:
        comm_out_path = out_dir + run_name + "_" + dna_int_graph.nodes[target_node_name]["in_name"].split(":")[0][0:80] + "_Top"
    cur_neighborhoods = community.partition_at_level(dendogram_com, level)
    cur_nodes = [nodes for nodes in cur_neighborhoods.keys() if cur_neighborhoods[nodes] == cur_neighborhoods[target_node_name]]
    source_community_graph = dna_int_graph.subgraph(cur_nodes)
    target_name = dna_int_graph.nodes[target_node_name]["in_name"].split(":")[0]
    source_community_graph.graph["target_gene"] = target_name
    
    orig_edges = source_community_graph.edges()
    if add_open:
        print("Annotating Open Regions")

        source_community_graph = add_openRegions_predict(source_community_graph, homo_gen, open_peaks_file, target_node_name, level)
        sys.stdout.flush()
        if add_snp:
            print("Annotating SNPs/CNVs")
            
            source_community_graph = add_variants_predict(source_community_graph, homo_gen, picked_tissue, target_node_name,out_dir)
            sys.stdout.flush()
    
    #Color Insulated Neighborhood according to node attributes, make mid distance to target node add min/max
    cur_chr = ""
    min_comm = 1e10
    max_comm = 0 
    for node in source_community_graph.nodes(data=True):
        #Keep min/max
        node_arr = re.split(r"[-:]",node[0])
        cur_chr = node_arr[0]
        #Set up colors
        source_community_graph.nodes[node[0]]["style"] = "filled"
        if (node[0] in target_nodes):
            source_community_graph.nodes[node[0]]["color"] = "rgb(17,166,216)"
            source_community_graph.nodes[node[0]]["size"] = 50
        elif ("gene" in node[1]):
            source_community_graph.nodes[node[0]]["color"] = "rgb(13,59,224)"
            source_community_graph.nodes[node[0]]["size"] = 25
        elif node[0] in enh_nodes:
            source_community_graph.nodes[node[0]]["color"] = "rgb(242,12,50)"
            source_community_graph.nodes[node[0]]["size"] = 25
        elif node[0] in ctcf_nodes:
            source_community_graph.nodes[node[0]]["color"] = "rgb(6,188,3)"
            source_community_graph.nodes[node[0]]["size"] = 25
        elif "top_open_snp" in node[1]:
            source_community_graph.nodes[node[0]]["color"] = "rgb(255,0,233)"
            source_community_graph.nodes[node[0]]["size"] = 12
        elif "dist2anchor" in node[1]:
            if (min(int(node_arr[1]),int(node_arr[2])) <= min_comm):
                min_comm = min(int(node_arr[1]),int(node_arr[2]))
            if (max(int(node_arr[1]),int(node_arr[2])) >= max_comm):
                max_comm = max(int(node_arr[1]),int(node_arr[2]))
            source_community_graph.nodes[node[0]]["color"] = "rgb(170,15,226)"
            source_community_graph.nodes[node[0]]["size"] = 15
        elif "all_ref" in node[1]:
            source_community_graph.nodes[node[0]]["color"] = "rgb(232,120,9)"
            source_community_graph.nodes[node[0]]["size"] = 10
        else:
            source_community_graph.nodes[node[0]]["color"] = "rgb(174,183,180)"
            source_community_graph.nodes[node[0]]["size"] = 20
    
    #Set min/max of nodes in community
    source_community_graph.graph["community_chr"] = cur_chr
    source_community_graph.graph["community_min"] = min_comm
    source_community_graph.graph["community_max"] = max_comm
    
    #get x,y,z local scatter_graphcoordinates and add to node
    position=nx.spring_layout(source_community_graph,dim=3)    
    for node in position:
        source_community_graph.nodes[node]["x"] = str(round(position[node][0],5))
        source_community_graph.nodes[node]["y"] = str(round(position[node][1],5))
        source_community_graph.nodes[node]["z"] = str(round(position[node][2],5))
    #Add UCSC session
#    source_community_graph.graph["ucsc_session"] = ucsc_session_in

    #If there's a previous graph, load and add it to current
    if os.path.exists(comm_out_path + "_graph.json"):
        with open(comm_out_path + "_graph.json","r") as json_comm:
            json_in = json.load(json_comm)
        source_community_graph_pre = json_graph.node_link_graph(json_in)
        source_community_graph = nx.compose(source_community_graph_pre,source_community_graph)
    #JSON output
    json_out_graph = json_graph.node_link_data(source_community_graph)
    json_comm = open(comm_out_path + "_graph.json","w+")
    json_out = json.dumps(json_out_graph)
    json_comm.write(json_out)
    json_comm.close()
    
    #Render HTML plotly output if lvl 0
    if level == 0:
        in_json_out_plotly_fig(comm_out_path + "_graph.json", 
                            make_3D=True, global_in=False, html_out=comm_out_path + "_3D_graph.html")
        in_json_out_plotly_fig(comm_out_path + "_graph.json", 
                            make_3D=False, global_in=False, html_out=comm_out_path + "_2D_graph.html")
    
    #2D Genome browser output
    for node in source_community_graph.nodes(data=True):
        if add_open and "dist2anchor" in node[1]:
            comm_open_out.write("\t".join(re.split(r"[-:]",node[0])) + "\n")
        elif add_snp and "all_ref" in node[1]:
            var_out = node[1]["gen_chrom"] + "\t" + str(node[1]["gen_start"]-1) + "\t" + str(node[1]["gen_end"]) + "\t" + node[1]["all_ref"] + "-->" + node[1]["all_alt"] 
            if "deep_score" in node[1]:
                var_out += ":" + str(node[1]["deep_score"]).replace(" ", "")
            comm_var_out.write(var_out + "\n")
        else:
            comm_nodes_out.write("\t".join(re.split(r"[-:]",node[0])) + "\n")
    for edge in orig_edges:
        if (edge2_strBed(edge)):
            comm_edges_out.write(edge2_strWashu(edge) + "\n")
            comm_edges_out_bed.write(edge2_strBed(edge) + "\n")
    comm_edges_out.write(edge2_strWashu(ctcf2bounds[target_node_name]) + "\n")
    comm_edges_out_bed.write(edge2_strBed(ctcf2bounds[target_node_name]) + "\n")
    return source_community_graph

def add_openRegions_predict(G, homo_gen, open_peaks, target_node_name, in_level):
    #Function to add open chromatin accessibility peaks to graph and predict TF binding using trained ML model
    #Making a new copy of the graph in order to add nodes
    G = G.copy()
    #how far can the chromatin accessibility be from the node to be considered?
    open_region_slop = 50000
    orig_edges = list(G.edges())
    orig_nodes = list(G.nodes())
    with open(open_peaks,"r") as all_peaks:
        for peak in all_peaks:
            arr = peak.strip().split("|")
            if len(arr) > 1 and (arr[1] != "NA") and (abs(int(arr[2])) <= open_region_slop):
                dist2anchor = abs(int(arr[2]))
                arr_atac = arr[0].strip().split()
                node_arr = arr[1].strip().split()
                atac_node = arr_atac[0] + ":" + str(int(arr_atac[1])) + "-" + str(int(arr_atac[2]))
                #Spoof open_size of peak (replace for p-value or something)
                open_size = 1
#                open_size = float(arr_atac[3])
                anchor_node = node_arr[0] + ":" + str(int(node_arr[1])) + "-" + str(int(node_arr[2]))
                if anchor_node in orig_nodes:
                    if in_level < 0:
                        deep_open_pred_out = deepbind_predict_tf_range(arr_atac[0],int(arr_atac[1]),int(arr_atac[2]),homo_gen, TF_RBP_ids, target_node_name)
                        #bass_open_pred_out = basset_predict_tf_range(arr_atac, homo_gen, target_node_name, openpeak2bassetpeak)
                    else:
                        deep_open_pred_out = ["","",""]
                    top_tf_names = deep_open_pred_out[0]
                    top_tf_probs = deep_open_pred_out[1]
                    top_tf_probs = [str(round(prob,2)) for prob in top_tf_probs]
                    top_tf_locs = deep_open_pred_out[2]
                    predict_tf_deepbind = ""
                    # Show top 15 TFs predicted under each chromatin accessibility peak
                    for tf in range(0,min([len(top_tf_names),15])):
                        predict_tf_deepbind += str(top_tf_names[tf]) + " [" + str(int(top_tf_locs[tf])) + "] :" + top_tf_probs[tf] + ", "
                    atac_mid = (int(arr_atac[1]) + int(arr_atac[2]))/2.0
                    top_tf_names=("|".join(top_tf_names))
                    top_tf_probs=("|".join(top_tf_probs))
                    filt_tf_open = top_tf_names[0:min([len(top_tf_names),15])]
                    G.add_node(atac_node, mid=atac_mid,name=atac_node, dist2anchor=dist2anchor, deepbind_tf=predict_tf_deepbind,top_tf_names=top_tf_names, filt_tf_open = filt_tf_open, top_tf_probs=top_tf_probs, open_size=open_size)
                    G.add_edge(anchor_node,atac_node,label="",weight=1.0)
    return G

def add_variants_predict(G, homo_gen, tissue_type, target_node_name,out_dir):
    #Function to add variants and predict their effect using DeepBind
    in_genes = G.nodes[target_node_name]["in_name"].split(",")
    #get target gene for eQTL analysis
    mg = mygene.MyGeneInfo()
    target_gene_out = mg.getgenes(in_genes, fields="ensembl", species="human", fetch_all=True)
    ensembl_target = ""
    for gene in target_gene_out:
        if "ensembl" in gene:
            ensembl_target = gene["ensembl"]["gene"]
    # flag for if any SNPs where added to the Graph
    add_any_SNPs = False
    #Save vcf output
#    basset_name_out = os.getcwd() + "/" + target_node_name.replace(":","_").replace("-","_") + "_basset"
#    if not os.path.exists(basset_name_out):
#        os.makedirs(basset_name_out)
#    open_vcf = open(basset_name_out + "/open_variants.vcf","w+")
    #find all variants that fall within these open peaks using MyVariant.info 
    mv = myvariant.MyVariantInfo()
    
    for node in list(G.nodes(data=True)):
        if "dist2anchor" in node[1]:
            #Pull all variants within chromatin accessibility peak fetch_all=True, 
            all_snps = mv.query(node[0], fetch_all=True, verbose=False)
            #limit number of SNPs under each peak
            snp_cnt = 0
            snp_limit = 10
            for snp in all_snps:
                if "dbsnp" not in snp:
                    continue
                
                if (snp_cnt >= snp_limit):
                    continue
                snp_cnt += 1
                
                add_any_SNPs = True      
                #Add CNVs/SNVs to graph
                G.add_node(snp["_id"],gen_start=snp["hg19"]["start"], all_ref = str(snp["dbsnp"]["ref"]),
                           all_alt = str(snp["dbsnp"]["alt"]),gen_end=snp["hg19"]["end"], gen_chrom="chr"+str(snp["chrom"]), 
                           name=snp["_id"])
                G.add_edge(snp["_id"], node[0],label="",weight=1.0)
                #print to VCF
#                open_vcf.write("chr"+str(snp["chrom"]) + "\t" + str(snp["hg19"]["start"]) + "\t" + snp["_id"] +
#                               "\t" + str(snp["dbsnp"]["ref"]) + "\t" + str(snp["dbsnp"]["alt"]) + "\n")
                
                snp_loc = str(G.nodes[snp["_id"]]["gen_chrom"]) + ":" + str(G.nodes[snp["_id"]]["gen_start"]) + "-" + str(G.nodes[snp["_id"]]["gen_end"])

#                    G.node[snp["_id"]]["DMC_node"] = DMC_var
                if ("cadd" in snp):
                    G.nodes[snp["_id"]]["cadd"] = snp["cadd"]["rawscore"]
                #add rsid if present and add GTEX eQTL pval and if target gene has Ensembl ID
                if (ensembl_target is not None) and ("dbsnp" in snp) and ("rsid" in snp["dbsnp"]):
                    serv = "http://rest.ensembl.org/eqtl/variant_name/homo_sapiens/"
                    ext = snp["dbsnp"]["rsid"].strip() + "?content-type=application/json;statistic=p-value;stable_id=" + ensembl_target + ";tissue=" + tissue_type 
                    r = requests.get(serv + ext, headers={ "Content-Type" : "application/json"})
                    if r.ok:
                        ret_val = r.json()
                        if isinstance(ret_val,list) and (len(ret_val)>0):
                            eqtl_pval = ret_val[0]["value"]
                            G.nodes[snp["_id"]]["gtex_eqtl_pval"] = round(eqtl_pval,6)
                    G.nodes[snp["_id"]]["rsid"] = snp["dbsnp"]["rsid"]
                #add GRASP phenotype and PMID information for SNP
                if ("grasp" in snp) and ("publication" in snp["grasp"]):
                    pheno_gwas = ""
                    pmid_gwas = ""
                    #single or multiple publications w/ SNP?
                    if isinstance(snp["grasp"]["publication"], list):
                        for pub in snp["grasp"]["publication"]:
                            if isinstance(pub["paper_phenotype_description"],list):
                                pheno_gwas += ",".join([str(pheno) + ", " for pheno in pub["paper_phenotype_description"]]) + ", "
                            else:
                                pheno_gwas += pub["paper_phenotype_description"] + ","
                            if isinstance(pub["pmid"],list):
                                pmid_gwas += ",".join([str(pmid) for pmid in pub["pmid"]]) + ", "
                            else:
                                pmid_gwas += str(pub["pmid"]) + ","
                    else:
                        pub = snp["grasp"]["publication"]
                        if isinstance(pub["paper_phenotype_description"],list):
                            pheno_gwas += "".join([str(pheno) + ", " for pheno in pub["paper_phenotype_description"]]) + ", "
                        else:
                            pheno_gwas += pub["paper_phenotype_description"]
                        if isinstance(pub["pmid"],list):
                            pmid_gwas += ",".join([str(pmid) for pmid in pub["pmid"]])
                        else:
                            pmid_gwas += str(pub["pmid"])
                    G.nodes[snp["_id"]]["grasp_pheno"] = pheno_gwas.replace(" ","_")
                    G.nodes[snp["_id"]]["grasp_pmid"] = pmid_gwas.replace(" ","_")
#    open_vcf.close()
    #predict effect of all variation added to graph using DeepBind and Basset
#    if add_any_SNPs:
#        print("Deepbind SNP prediction of: " + target_node_name)
#        G = deepbind_predict_SNPs(G, homo_gen, TF_RBP_ids,target_node_name)
#        print("Basset SNP prediction of: " + target_node_name)
#        G = add_basset_sad_sat(G, homo_gen, target_node_name,out_dir)
#    try:
#        os.remove(basset_name_out + "open_variants.vcf")
#    except OSError:
#        pass
#    shutil.rmtree(basset_name_out)
    return G

def build_IN_levels(target_node_name, dendogram_com, dna_int_graph, genome_fa, out_dir, run_name, json_draft_name, open_peaks_file,picked_tissue, ctcf2bounds,rough_3D_loc):
    
    #Genome browser output tracks
    #Edges in washu format
    comm_edges_out = open(out_dir + run_name + "_target_community_edges_lvl_1.washu","a+")
    comm_edges_out_lvl1 = open(out_dir + run_name + "_target_community_edges_lvl_2.washu","a+")
    comm_edges_out_lvl2 = open(out_dir + run_name + "_target_community_edges_lvl_Top.washu","a+")
    #Edges in bed format
    comm_edges_out_bed = open(out_dir + run_name + "_target_community_edges_lvl_1.bed","a+")
    comm_edges_out_lvl1_bed = open(out_dir + run_name + "_target_community_edges_lvl_2.bed","a+")
    comm_edges_out_lvl2_bed = open(out_dir + run_name + "_target_community_edges_lvl_Top.bed","a+")
    #Nodes in bed format
    comm_nodes_out = open(out_dir + run_name + "_target_community_nodes_lvl_1.bed","a+")
    comm_nodes_out_lvl1 = open(out_dir + run_name + "_target_community_nodes_lvl_2.bed","a+")
    comm_nodes_out_lvl2 = open(out_dir + run_name + "_target_community_nodes_lvl_Top.bed","a+")
    #Variants in bed format
    comm_var_out = open(out_dir + run_name + "_target_community_variants_lvl_1.bed","a+")
    comm_var_out_lvl1 = open(out_dir + run_name + "_target_community_variants_lvl_2.bed","a+")
    comm_var_out_lvl2 = open(out_dir + run_name + "_target_community_variants_lvl_Top.bed","a+")
    #Open regions in bed format
    comm_open_out = open(out_dir + run_name + "_target_community_open_regions_lvl_1.bed","a+")
    comm_open_out_lvl1 = open(out_dir + run_name + "_target_community_open_regions_lvl_2.bed","a+")
    comm_open_out_lvl2 = open(out_dir + run_name + "_target_community_open_regions_lvl_Top.bed","a+")


    target_name = dna_int_graph.nodes[target_node_name]["in_name"].split(":")[0]
    chr_in = re.split(r"[-:]",target_node_name)[0]
#    ucsc_session_in = ucsc_session + "&position=" + chr_in + ":" + str(dna_int_graph.node[target_node_name]["in_max"]) + "-" + str(dna_int_graph.node[target_node_name]["in_min"])
    
    print("Annotating Insulated Neighborhood: " + target_name)
    #load genome from genome fasta
    homo_gen = pysam.FastaFile(genome_fa)
    #fetch dendogram at CTCF-CTCF loop restricted resolution
    dendogram_com = target2partition[target_node_name]
    
    #level 0
    level = 0
    add_open = True
    if dna_int_graph.graph["species"] == "Human":
        add_snp = True
    else:
        add_snp = False
    source_community_graph = annotate_IN(target_node_name, dendogram_com, level, dna_int_graph, add_open, add_snp,
                                                  homo_gen, open_peaks_file, picked_tissue, ctcf2bounds, comm_open_out,
                                                  comm_var_out, comm_nodes_out, comm_edges_out, comm_edges_out_bed, out_dir)
    
    #level 1 
    if len(dendogram_com) > 1:
        level = 1
        add_open = False
        add_snp = False
        source_community_graph_lvl1 = annotate_IN(target_node_name, dendogram_com, level, dna_int_graph, add_open, add_snp,
                                                  homo_gen, open_peaks_file, picked_tissue, ctcf2bounds, comm_open_out_lvl1,
                                                  comm_var_out_lvl1, comm_nodes_out_lvl1, comm_edges_out_lvl1, comm_edges_out_lvl1_bed, out_dir)
        
    #top level
    if len(dendogram_com) > 2:
        level = len(dendogram_com)-1
        add_open = False
        add_snp = False
        source_community_graph_lvl2 = annotate_IN(target_node_name, dendogram_com, level, dna_int_graph, add_open, add_snp,
                                                  homo_gen, open_peaks_file, picked_tissue, ctcf2bounds, comm_open_out_lvl2,
                                                  comm_var_out_lvl2, comm_nodes_out_lvl2, comm_edges_out_lvl2, comm_edges_out_lvl2_bed, out_dir)
        
    #Add rough XYZ locations for neighborhoods, not for ChrX,Y,M,22   
    comm_out_path = run_name + "_" + dna_int_graph.nodes[target_node_name]["in_name"].split(":")[0][0:80] + "_lvl1"
    if "chrX" not in target_node_name and "chrY" not in target_node_name and "chr22" not in target_node_name and "chrM" not in target_node_name:
        with open(json_draft_name,"r") as json_comm:
            json_in = json.load(json_comm)
        draft_genome_out = json_graph.node_link_graph(json_in)
        draft_genome_out.add_node(target_node_name, name=target_name, json_name=(comm_out_path + "_graph.json"), edges_ucsc_lvl0 = run_name + "_target_community_edges_lvl_1.bed", nodes_lvl0_bed = run_name+"_target_commmunity_nodes_lvl1.bed",
	 open_lvl0_bed = run_name + "_target_community_open_lvl_1.bed", var_lvl0_bed = run_name + "_target_community_variants_lvl_1.bed")
        #add HSA derived XYZ coordinates for IN collections
        draft_genome_out = add_xyz_loc(draft_genome_out,rough_3D_loc)
        json_out_graph = json_graph.node_link_data(draft_genome_out)
        json_out = json.dumps(json_out_graph)
        with open(json_draft_name,"w") as json_comm:
            json_comm.write(json_out)
        
    comm_open_out.close()
    comm_var_out.close()
    comm_nodes_out.close()
    comm_edges_out.close()
    comm_edges_out_bed.close()
    comm_open_out_lvl1.close()
    comm_var_out_lvl1.close()
    comm_nodes_out_lvl1.close()
    comm_edges_out_lvl1.close()
    comm_edges_out_lvl1_bed.close()
    comm_open_out_lvl2.close()
    comm_var_out_lvl2.close()
    comm_nodes_out_lvl2.close()
    comm_edges_out_lvl2.close()
    comm_edges_out_lvl2_bed.close()
    
def add_xyz_loc(G, xyz_bed):
    sub_nodes = str(uuid.uuid4()) + ".txt"
    G_nodes = open(sub_nodes,"w+")
    for node in G.nodes():
        arr = re.split(r"[-:]",node)
        G_nodes.write(arr[0] + "\t" + str(arr[1]) + "\t" + str(arr[2]) + "\n")
    G_nodes.close()
    xyz_nodes = str(uuid.uuid4()) + ".txt"
    os.system("sort-bed " + sub_nodes + " | closest-features --closest - "+ xyz_bed + " > " + xyz_nodes)
    with open(xyz_nodes,"r") as xyz_nodes_file:
        for global_coord_node in xyz_nodes_file:
            arr = global_coord_node.strip().split("|")
            graph_node = arr[0].split("\t")
            xyz_node = arr[1].split("\t")
            id_node = graph_node[0]+":"+graph_node[1]+"-"+graph_node[2]
            global_pos = (xyz_node[3]).split(",")
            G.nodes[id_node]["x"] = global_pos[0]
            G.nodes[id_node]["y"] = global_pos[1]
            G.nodes[id_node]["z"] = global_pos[2]
    os.remove(sub_nodes)
    os.remove(xyz_nodes)
    return G
             
        
#plotly plotting function for networkx graphs
def in_json_out_plotly_fig(json_name, make_3D, global_in, html_out):
    #Load serialized JSON and load into networkx graph
    with open(json_name,"r") as json_comm:
        json_in = json.load(json_comm)
        G = json_graph.node_link_graph(json_in)
    
    #Set up axis features for plotting
    axis=dict(
              showline=False,
              zeroline=False,
              showgrid=False,
              showticklabels=False,
              showbackground=False,
              ticks='',
              title="")
    #Set up layout for plotly plotting
    if global_in:
        layout = go.Layout( scene = dict(
                                         xaxis=axis,
                                         yaxis=axis,
                                         zaxis=axis),
                           plot_bgcolor='rgba(0,0,0,0)',
                           width=720,
                           height=720,
                           showlegend=False,
                           margin=dict(l=0,
                                       r=0,
                                       b=0,
                                       t=20),
                           hovermode="closest",
                           title=go.layout.Title(text="Global View", x=.5, y=.97)
                          )
    else:
        layout = go.Layout( scene = dict(
                                         xaxis=axis,
                                         yaxis=axis,
                                         zaxis=axis),
                            plot_bgcolor='rgba(0,0,0,0)',
                            showlegend=False,
                            width=1280,
                            height=720,
                            margin=dict(l=0,r=0,b=0,t=20),
                            hovermode="closest",
                            title=go.layout.Title(text="Target neighborhood: " + G.graph["target_gene"],
                                                 x=.5, y=.97)
                          )
        
    #make node trace
    node_text = []
    node_color = []
    node_size = []
    node_opacity = []
    node_x = []
    node_y = []
    #Make 3D?
    if make_3D:
        node_z = []
    for node in G.nodes(data=True):
        text_node = ""
        if "rsid" in node[1]:
            text_node = "rsid: " + node[1]["rsid"] + "\n" + node[1]["name"]
        elif "name" in node[1] and "fpkm" in node[1]:
            text_node += node[0] + " " + node[1]["name"]
        elif "name" in node[1]:
            text_node += node[1]["name"]
        else:
            text_node += node[0]
        if "fpkm" in node[1]:
            text_node += "\nFPKM: " + str(node[1]["fpkm"])
        if "tfbs" in node[1]:
            text_node += "\n" + textwrap.fill("TFBS: " + node[1]["tfbs"])
        if "cadd" in node[1]:
            text_node += "\n" + textwrap.fill("CADD Score: " + str(node[1]["cadd"]))
        if "deep_score" in node[1]:
            text_node += "\n" + textwrap.fill("DeepBind Delta: " + str(node[1]["deep_score"]))
        if "deepbind_tf" in node[1]:
            text_node += "\n" + textwrap.fill("DeepBind Top TFs: " + str(node[1]["deepbind_tf"]))
        if "grasp_pheno" in node[1]:
            text_node += "\n" + textwrap.fill("GWAS Phenotype (GRASP): " + str(node[1]["grasp_pheno"]))
        if "grasp_pmid" in node[1]:
            text_node += "\n" + textwrap.fill("GWAS PMID (GRASP): " + str(node[1]["grasp_pmid"]))
        if "gtex_eqtl_pval" in node[1]:
            text_node += "\n" + textwrap.fill("GTEX eQTL P-value: " + str(node[1]["gtex_eqtl_pval"]))
        if "sad_abs_sum" in node[1]:
            text_node += "\n" + textwrap.fill("Sum of SAD across tissues: " + str(node[1]["sad_abs_sum"]))
        node_text.append(text_node.replace("\n","<br>"))
        node_x.append(node[1]["x"])
        node_y.append(node[1]["y"])
        if make_3D:
            node_z.append(node[1]["z"])
        if "color" in node[1]:
            node_color.append(node[1]["color"])
        else:
            node_color.append("white")
        if "size" in node[1]:
            node_size.append(node[1]["size"])
        else:
            node_size.append(10)
        if "opacity" in node[1]:
            node_opacity.append(node[1]["opacity"])
        else:
            node_opacity.append(1)

    if make_3D:
        traceN = go.Scatter3d(x=node_x, 
                              y=node_y, 
                              z=node_z, 
                              mode="markers", 
                              text=node_text,
                              marker=go.scatter3d.Marker(color=node_color,
                                                       size=node_size,
                                                       opacity=1), 
                              name='', 
                              hoverinfo='text')
    else:
        traceN = go.Scatter(x=node_x, 
                    y=node_y, 
                    mode="markers", 
                    text=node_text, 
                    marker=go.scatter.Marker(color=node_color,
                                             size=node_size,
                                             opacity=node_opacity), 
                    name='', 
                    hoverinfo='text')
    #make edge trace and annotations
    edges_x = []
    edges_y = []
    edges_txt = []
    edges_annot_x = []
    edges_annot_y = []
    edges_annot_color = []
    edges_annot_size = []
    if make_3D:
        edges_z = []
        edges_annot_z = []
        
    for edge in G.edges(data=True):
        edges_x += [G.nodes[edge[0]]["x"],G.nodes[edge[1]]["x"], None]
        edges_y += [G.nodes[edge[0]]["y"],G.nodes[edge[1]]["y"], None]
        if make_3D:
            edges_z += [G.nodes[edge[0]]["z"],G.nodes[edge[1]]["z"], None]
        text_edge = ""
        if "weight" in G.get_edge_data(*edge):
            text_edge += "Weight: " + str(G.get_edge_data(*edge)["weight"])
            edges_annot_x.append((float(G.nodes[edge[0]]["x"])+float(G.nodes[edge[1]]["x"]))/2.0)
            edges_annot_y.append((float(G.nodes[edge[0]]["y"])+float(G.nodes[edge[1]]["y"]))/2.0)
            if make_3D:
                edges_annot_z.append((float(G.nodes[edge[0]]["z"])+float(G.nodes[edge[1]]["z"]))/2.0)               
            edges_annot_color.append("black")
            edges_annot_size.append(5)
        if text_edge:
            edges_txt.append(text_edge)
        
    if make_3D:
        traceE = go.Scatter3d(x=edges_x, 
                              y=edges_y, 
                              z=edges_z, 
                              mode="lines", 
                              hoverinfo = "none")
        traceE_annot = go.Scatter3d(x=edges_annot_x, 
                                    y=edges_annot_y, 
                                    z=edges_annot_z, 
                                    name="", 
                                    mode="markers", 
                                    text=edges_txt,
                                    hoverinfo="text",
                                    marker=go.scatter3d.Marker(color=edges_annot_color,
                                                     size=edges_annot_size,
                                                     opacity=.5),
                                   )
    else:
        traceE = go.Scatter(x=edges_x, 
                    y=edges_y, 
                    name="", 
                    mode="lines", 
                    hoverinfo = "none")
        traceE_annot = go.Scatter(x=edges_annot_x, 
                                  y=edges_annot_y, 
                                  name="", 
                                  mode="markers", 
                                  text=edges_txt,
                                  hoverinfo="text",
                                  marker=go.scatter.Marker(color=edges_annot_color,
                                                           size=edges_annot_size,
                                                           opacity=.5),
                                 )
        
    #Make plotly figure object
    data = [traceE, traceN, traceE_annot]
    fig = go.Figure(data=data, layout=layout)

    #Remove tick labels because for whatever reason not working in layout dictionary
    fig.update_xaxes(showticklabels=False)
    fig.update_yaxes(showticklabels=False)
    
    #Write rendered plotly figure to html output
    pio.write_html(fig, file=html_out, auto_open=False)

#2D plotting function using pygenometracks take in annotated global JSON generate all images
def in_json_out_genome_views(global_json_name, view_out, track_dir, IN_dir, run_name):
    #Load serialized JSON and load into networkx graph
    with open(IN_dir+global_json_name,"r") as json_comm:
        #json load is incorrectly adding a bunch of nodes, no idea why
        json_in = json.load(json_comm)
        G = json_graph.node_link_graph(json_in)
    
    #Get needed features of graph
    cell_type = G.graph["cell_type"]
    TAD_regions = G.graph["TAD_regions"]
    tracks_ini_fine = G.graph["tracks_ini_fine"]
    tracks_ini_broad = G.graph["tracks_ini_broad"]    
    #Make 4 zooms and output to file for viewing        
    # 4 zoom levels: TSS, gene, IN, TAD
    tss_zoom = view_out + run_name + "_tss_zoom"
    tss_zoom_in_bed = open(tss_zoom + ".bed", "w")
    gene_zoom = view_out + run_name + "_gene_zoom"
    gene_zoom_in_bed = open(gene_zoom + ".bed", "w")
    IN_zoom = view_out + run_name + "_IN_zoom"
    IN_zoom_in_bed = open(IN_zoom + ".bed", "w")
    TAD_zoom = view_out + run_name + "_TAD_zoom"
    TAD_zoom_in_bed = open(TAD_zoom + ".bed", "w")

    for IN_node in G.nodes(data=True):
        #Get IN name from graph
        in_name = IN_node[1]["name"].replace(" ",":")
        
        #Pull in target IN json and get target TSS node
        with open(IN_dir+IN_node[1]["json_name"],"r") as json_comm:
            json_in = json.load(json_comm)
            G_in = json_graph.node_link_graph(json_in)
        cur_target = G_in.nodes[IN_node[0]]
        
        #Fetch gene annotation from mygene.info
        mg = mygene.MyGeneInfo()
        target_gene_out = mg.querymany(cur_target["name"].split(","),
                                       scopes=["ensemblgene", "symbol"],
                                   fields=["genomic_pos"], verbose=False,
                                   species="human", fetch_all=True)
        
        #Get TSS limits from target node
        gene_arr = re.split(r"[-:]",IN_node[0])
        #TSS limits
        cur_cur = gene_arr[0]
        tss_min = gene_arr[1]
        tss_max = gene_arr[2]
        tss_zoom_in_bed.write(cur_cur + "\t" +
                              str(tss_min) + "\t" +
                              str(tss_max) + "\t" +
                              in_name + "\n")
        
        #Gene zoom limits
        gene_min = 1e10
        gene_max = 0 
        for gene in (target_gene_out):
            if (("genomic_pos" in gene) 
                and ("start" in gene["genomic_pos"]) 
                and ("end" in gene["genomic_pos"])):
                if (max(int(gene["genomic_pos"]["start"]),
                       int(gene["genomic_pos"]["end"])) > gene_max):
                    gene_max = max(int(gene["genomic_pos"]["start"]),
                               int(gene["genomic_pos"]["end"]))
                if (min(int(gene["genomic_pos"]["start"]),
                       int(gene["genomic_pos"]["end"])) < gene_min):
                    gene_min = min(int(gene["genomic_pos"]["start"]),
                               int(gene["genomic_pos"]["end"]))
        gene_zoom_in_bed.write(cur_cur + "\t" +
                              str(gene_min) + "\t" +
                              str(gene_max) + "\t" +
                              in_name + "\n")

        #IN limits
        cur_min = G_in.graph["community_min"]
        cur_max = G_in.graph["community_max"]
        IN_zoom_in_bed.write(cur_cur + "\t" +
                              str(cur_min) + "\t" +
                              str(cur_max) + "\t" +
                              in_name + "\n")        

    #Close and sort TSS file
    tss_zoom_in_bed.close()
    os.system("sort-bed " + tss_zoom + ".bed > t.bed; mv t.bed " + tss_zoom + ".bed")

    #Get closest TAD to each TSS
    TAD_slop = 50000
    tad_command = ("closest-features --dist --closest " + tss_zoom + ".bed " +
              TAD_regions + " > TAD_closest.bed")
    os.system(tad_command)
    #Assign each TSS to closest TAD
    with open("TAD_closest.bed") as closest_TADs:
        for TAD in closest_TADs:
            arr = TAD.strip().split("|")
            if (len(arr)>1) and (arr[1] != "NA") and (abs(int(arr[2])) <= TAD_slop) :
                tss_arr = arr[0].split()
                tad_arr = arr[1].split()
                TAD_zoom_in_bed.write(str(tad_arr[0]) + "\t" +
                                      str(tad_arr[1]) + "\t" +
                                      str(tad_arr[2]) + "\t" +
                                      in_name + "\n")  
    
    #Run 4 runs of pygenometracks for each zoom level
    tss_zoom_in_bed.close()
    gene_zoom_in_bed.close()
    IN_zoom_in_bed.close()
    TAD_zoom_in_bed.close()
    
    #Collapse TSS,gene into fine
    os.system("cat " + tss_zoom + ".bed " + gene_zoom + ".bed | sort-bed - > fine_in.bed")
    viz_command = ("pyGenomeTracks --tracks " + tracks_ini_fine + 
                    " --BED fine_in.bed --dpi 500 " +
                    " --outFileName " + view_out + run_name + "_fine_zoom.png")
    print("Generating fine zoom genomic views")    
    os.system(viz_command)
    
    #Collapse IN, TAD into broad
    os.system("cat " + IN_zoom + ".bed " + TAD_zoom + ".bed | sort-bed - > broad_in.bed")
    viz_command = ("pyGenomeTracks --tracks " + tracks_ini_broad + 
                    " --BED broad_in.bed --dpi 500 " +
                    " --outFileName " + view_out + run_name + "_broad_zoom.png")
    print("Generating broad zoom genomic views")    
    os.system(viz_command)
    #os.remove("TAD_closest.bed")
    os.remove("broad_in.bed")
    os.remove("fine_in.bed")               
    
if __name__ == "__main__":
    #Handle command line arguments
    parser = argparse.ArgumentParser(description="Pipeline to output graphical genomic model from 2D data")
    parser.add_argument("-i","--input_directory", help="Where all the input data lives", required=True)
    parser.add_argument("-a","--annotation_directory", help="Where all the annotations (ex: hg19.fa) live", required=True)
    parser.add_argument("-o","--output_directory", help="Output directory of graphs", required=True)
    parser.add_argument("-t","--insulated_neigh", help="Comma separated list of Insulated Neighborhoods to build", required=False)
    parser.add_argument("-c","--dna_int_cell_type", help="Cell type to use in building Insulated Neighborhoods", required=True)
    parser.add_argument("-q","--eQTL_tissue", help="What tissue from GTeX to do eQTL analysis w/. Available tissues: Adipose_Subcutaneous," + 
                        "Adipose_Visceral_Omentum, Adrenal_Gland, Artery_Aorta, Artery_Coronary, Artery_Tibial, Brain_Anterior_cingulate_cortex_BA24," +
                        " Brain_Caudate_basal_ganglia, Brain_Cerebellar_Hemisphere, Brain_Cerebellum, Brain_Cortex, Brain_Frontal_Cortex_BA9, " +
                        "Brain_Hippocampus, Brain_Hypothalamus, Brain_Nucleus_accumbens_basal_ganglia, Brain_Putamen_basal_ganglia," +
                        "Breast_Mammary_Tissue, Cells_EBV-transformed_lymphocytes, Cells_Transformed_fibroblasts, Colon_Sigmoid, " +
                        "Colon_Transverse, Esophagus_Gastroesophageal_Junction, Esophagus_Mucosa, Esophagus_Muscularis, Heart_Atrial_Appendage, " +
                        "Heart_Left_Ventricle, Liver, Lung, Muscle_Skeletal, Nerve_Tibial, Ovary, Pancreas, Pituitary, Prostate, " +
                        "Skin_Not_Sun_Exposed_Suprapubic, Skin_Sun_Exposed_Lower_leg, Small_Intestine_Terminal_Ileum, Spleen, Stomach, " +
                        "Testis, Thyroid, Uterus, Vagina, Whole_Blood", required=True)
    parser.add_argument("-g","--top_gene_fpkm", help="Instead of defining Insulated Neighborhoods take genes expressed above FPKM cutoff", required=False)
    parser.add_argument("-n","--name_output", help="Output prefix for all output files", required=True)
    parser.add_argument("-m","--marker_region", help="Build all Insulated Neighborhoods of genes which fall within this chromosomal region (ex: chr1:54996039-55996039)", required=False)
    args = vars(parser.parse_args())
    
    cell_type = (args["dna_int_cell_type"]).lower() #Primed, Naive, Neuron, NPC, GM12878, Jurkatt, K562, HeLa, mESC
    picked_tissue = args["eQTL_tissue"]
    home_dir = args["input_directory"]
    ann_dir = args["annotation_directory"]
    run_name = args["name_output"]
    out_dir = args["output_directory"] + "/" + cell_type + "/" + run_name + "/"
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    #switch to output directory for analysis
    os.chdir(out_dir)

    input_dir = dict()
    if cell_type == "lung":
        genome_version = "GRCh38"
        species = "Human"
        rough_3D_loc = home_dir + "hIMR90_miniMDS_hg38.bed"
        dna_ints = home_dir + "primed_.7_origami_intra_hg38.arcs" 
        enh_peaks = home_dir + "ENCFF731CHG_H3k27Ac_NHLF_hg38.bed" 
        prom_peaks = home_dir + "ENCFF916HZY_H3k4Me3_SAEC_hg38.bed" 
        silenced_peaks = home_dir + "ENCFF959LWN_H3k9Me3_NHLF_hg38.bed"
        TAD_regions = home_dir + "ENCFF307RGV_IMR90_TAD_hg38.bed"
        ctcf_peaks = home_dir + "ENCFF777ODE_CTCF_NHLF.bed"
        open_peaks = home_dir + "ENCFF827VFY_DNAse_SAEC.bed" 
        rna_seq = home_dir + "ENCFF679UBR_HBEC_rnaseq_transcript.tsv" 
        tss_ensembl = ann_dir + "biomart_hg38_TSS_gene.txt"
#        encode_bed = ann_dir + "hg19_encode_TFBSv2.bed"
        genome_fa =  ann_dir + "GRCh38_filt_dna_rm.fa"
#        TF_RBP_ids = home_dir + "TF_RBP_human.ids"
        tracks_ini_fine = home_dir + "regulatory_architecture_lung_fine.ini"
        tracks_ini_broad = home_dir + "regulatory_architecture_lung_broad.ini"
        add_snp = True
    else:
        raise NameError("Cant find that cell type, try: lung")

    #get list of included chromosomes from genome fasta
    homo_gen = pysam.FastaFile(genome_fa)
    #keep only canonical names chr1-22,X,Y,M
    chr_names = homo_gen.references
    chr_names = [("chr"+chro) for chro in chr_names if "_" not in chro]
    #initiate genome graph
    dna_int_graph = nx.Graph(style="filled")
    #Set genomewide attributes
    dna_int_graph.graph["species"] = species
    dna_int_graph.graph["genome_version"] = genome_version
    dna_int_graph.graph["tissue_type"] = picked_tissue
    #Load DNA interactions
    with open(dna_ints) as dna_ints_iter:
        for dna_int in dna_ints_iter:
            arr = dna_int.split()
            x = str(arr[0]) + ":" + str(arr[1]) + "-" + str(arr[2])
            mid_x = (int(arr[1])+int(arr[2]))/2.0
            y = str(arr[3]) + ":" + str(arr[4]) + "-" + str(arr[5])
            mid_y = (int(arr[4])+int(arr[5]))/2.0
            if (arr[0] in chr_names) and (arr[3] in chr_names):
                dna_int_graph.add_edge(x,y,label=1, capacity = 1, weight=float(arr[6]))
    #set up nodes bed files
    dna_int_graph = uniq_nodes(dna_int_graph)
    #get classes of nodes
    h3k4me3_nodes = overlap_graph_bed(dna_int_graph, prom_peaks)
    enh_nodes = overlap_graph_bed(dna_int_graph, enh_peaks)
    ctcf_nodes = overlap_graph_bed(dna_int_graph, ctcf_peaks)
    #Label corresponding nodes as Genes and add fpkm attribute
    dna_int_graph = label_gene_nodes(dna_int_graph, rna_seq, tss_ensembl)
    #add TFBS from Chip-Seq under nodes
    #dna_int_graph = add_TFBS_graph(dna_int_graph, encode_bed)
    
    #get target nodes and target node names
    target_nodes = get_target_nodes(dna_int_graph, args, cell_type)
    target_communities = set()
    for node in target_nodes:
        target_communities |=  set((dna_int_graph.nodes[node]["name"]).split(","))
    print("Number of Target Genes: " + str(len(target_communities)))
    print("Target Genes to build Insulated Neighborhoods: " + str(",".join(target_communities)))
    #Sort nodes by FPKM for louvain
    node_fpkm = dict()
    start_comm = dict()
    for node in dna_int_graph.nodes(data=True):
        if ("name" in node[1]):
            node_fpkm[node[0]] = float(node[1]["fpkm"])
        else:                             
            node_fpkm[node[0]] = 0.0
    sorted_nodes = sorted(node_fpkm.keys(), key=lambda k: node_fpkm[k], reverse=True)
    node_count = 0
    for node in sorted_nodes:
        start_comm[node] = node_count
        node_count += 1

    #Set up CTCF-CTCF loop bed
    ctcf_loops = open("ctcf_loops.bed","w")
    for edge in dna_int_graph.edges(ctcf_nodes):
        l_node = re.split(r"[-:]",edge[0])
        r_node = re.split(r"[-:]",edge[1])
        if l_node[0] == r_node[0]:
            bounds = [int(l_node[1]),int(l_node[2]),int(r_node[1]),int(r_node[2])]
            ctcf_loops.write(l_node[0] + "\t" + str(int(min(bounds)))+"\t"+str(int(max(bounds))+1) + 
                             "\t" + edge[0] + ";" + edge[1] + "\n")
    ctcf_loops.close()
    os.system("sort-bed ctcf_loops.bed > t.bed; mv t.bed ctcf_loops.bed")   

    #Checking resolution parameter
    print("Checking resolution parameter against CTCF-CTCF loop")
    size_sum = list()
    size_avg = list()
    res_run = list()
    min_comm_nodes = dict()
    max_comm_nodes = dict()
    perc_CC = list()
    target2partition = dict()
    target2maxres = dict()
    ctcf2bounds = dict()
    num_res_steps = 5
    #dict of lists w/ sizes of community of target gene at all resolutions
    comm_target_lengths = dict()
    for i in np.logspace(-5,-1,num_res_steps):
        comm_out_name = out_dir + "min_max_communities_res_genes" + str(i) + ".bed"
        cur_comm_limits = open(comm_out_name, "w")
        res_run.append(i)
        dendogram_com = community.generate_dendrogram(dna_int_graph, resolution=i, part_init=start_comm)
        dna_int_comm = community.partition_at_level(dendogram_com, 0)
        size2comm = Counter(dna_int_comm.values())
        size_comms = [size for size in size2comm.values()]
        size_sum.append(len(set(dna_int_comm.values())))
        size_avg.append(sum(size_comms)/(1.0*len(size_comms)))
        #assess quality of communities generated
        #iterate over dna_int_comm keys and add each start/end to set stored in dict of comm --> min and comm-->max
        inv_map = {}
        for k, v in dna_int_comm.items():
            inv_map[v] = inv_map.get(v, [])
            inv_map[v].append(k)
        for comm, nodes in inv_map.items():
            #only print communities w/ genes in them
            gene_nodes = [node for node in nodes if "gene" in dna_int_graph.nodes[node]]
            target_overlap = [node for node in target_nodes if node in nodes]
            #get max and min of neighborhoods and write to disk
            if (len(gene_nodes)>0) and len(target_overlap)>0:
                if target_overlap[0] not in comm_target_lengths:
                    comm_target_lengths[target_overlap[0]] = list()
                comm_target_lengths[target_overlap[0]].append(str(len(nodes)))
                ref_node_arr = re.split(r"[-:]",nodes[0])
                min_comm = 1e10
                max_comm = 0
                for node in nodes:
                    node_arr = re.split(r"[-:]",node)
                    if node_arr[0] == ref_node_arr[0]:
                        if (min(int(node_arr[1]),int(node_arr[2])) <= min_comm):
                            min_comm = min(int(node_arr[1]),int(node_arr[2]))
                        if (max(int(node_arr[1]),int(node_arr[2])) >= max_comm):
                            max_comm = max(int(node_arr[1]),int(node_arr[2]))
                in_name_all_genes = ",".join(set([dna_int_graph.nodes[target_gene]["name"] for target_gene in target_overlap]))
                for target_gene in target_overlap:
                    dna_int_graph.nodes[target_gene]["comm_num"] = comm
                    dna_int_graph.nodes[target_gene]["in_name"] = in_name_all_genes
                    dna_int_graph.nodes[target_gene]["in_min"] = min_comm
                    dna_int_graph.nodes[target_gene]["in_max"] = max_comm
                cur_comm_limits.write(ref_node_arr[0] + "\t" + str(min_comm) + "\t" + str(max_comm+1) + "\t" + target_overlap[0] + "\n")
        cur_comm_limits.close()
        os.system("sort-bed " + comm_out_name + " > t.bed; mv t.bed " + comm_out_name)
        #count how many of target communities fell within CTCF - CTCF loop
        os.system("bedops -e -1 ctcf_loops.bed " + comm_out_name + " | bedmap --echo --multidelim \"_\" --echo-map-id-uniq  --fraction-ref .8 " + comm_out_name + " - > target_comm_inCC.bed")
        contained = 0.0
        for line in open("target_comm_inCC.bed"):
            line_arr = line.strip().split("|")
            part_arr = line_arr[0].strip().split()
            if len(line_arr)>1 and line_arr[1] != "":
                loop_arr = line_arr[1].strip().split("_")
                weights_cc = np.argmax([dna_int_graph.get_edge_data(*(loop.split(";")))["weight"] for loop in loop_arr])
                #Only output strongest CTCF bounding loop
                ctcf2bounds[part_arr[3]] = [(loop.split(";")[0],loop.split(";")[1]) for loop in loop_arr][weights_cc]
                contained += 1.0
                target2partition[part_arr[3]] = dendogram_com
                target2maxres[part_arr[3]] = i
        perc_CC.append((contained/len(comm_target_lengths))*100.0)
        os.remove("target_comm_inCC.bed")
        os.remove(comm_out_name)
    
    #Graph which is collection of neighborhoods and make new file
    draft_genome_out = nx.Graph()
    
    #Set genome wide attributes
    draft_genome_out.graph["species"] = species
    draft_genome_out.graph["genome_version"] = genome_version
    draft_genome_out.graph["tissue_type"] = picked_tissue
    draft_genome_out.graph["cell_type"] = cell_type
    draft_genome_out.graph["TAD_regions"] = TAD_regions
    draft_genome_out.graph["tracks_ini_fine"] = tracks_ini_fine
    draft_genome_out.graph["tracks_ini_broad"] = tracks_ini_broad
    
    #Start global graph
    dna_int_graph.graph["top_genes"] = ""
    json_draft_name = (out_dir + run_name[0:80]+ "_global_graph.json")
    json_comm = open(json_draft_name,"w+")
    json_out_graph = json_graph.node_link_data(draft_genome_out)
    json_out = json.dumps(json_out_graph)
    json_comm.write(json_out)
    json_comm.close()

    #Get unique target list of communities being annotated
    final_in_names = set()
    final_in_nodes = set()
    final_comms = set()
    for target_node in target2partition.keys():
        in_name = dna_int_graph.nodes[target_node]["in_name"]
        comm_target = dna_int_graph.nodes[target_node]["comm_num"]
        if comm_target not in final_comms:
            final_comms.add(str(comm_target))
            final_in_names.add(in_name)
            final_in_nodes.add(target_node)
    final_in_names = set(final_in_names)
    final_in_nodes = set(final_in_nodes)
    print("Annotating " + str(len(final_in_names)) + " Neighborhoods")
    print("Annotating the following Insulated Neighborhoods: " + "|".join(final_in_names))
    print("Annotating the following communities: " + "|".join(final_comms))
    
    #Assign each open chromatin accessibility peak to the closest cohesin peak
    open_peaks_file="open_peaks_nodes.bed"
    os.system("closest-features --closest --dist " + open_peaks + " all_nodes_sorted.bed > " + open_peaks_file)
    
    #Annotate neighborhoods
    #Parallelized version
    num_cores = multiprocessing.cpu_count()-1
    num_cores = 4
    results = Parallel(n_jobs=num_cores,temp_folder=out_dir) (delayed(build_IN_levels) 
                                                              (target_node_name,
                                                               dendogram_com,
                                                               dna_int_graph,
                                                               genome_fa,
                                                               out_dir,
                                                               run_name,
                                                               json_draft_name,
                                                               open_peaks_file,
                                                               picked_tissue,
                                                               ctcf2bounds,
                                                               rough_3D_loc)
                                                              for target_node_name in final_in_nodes)
    
    #Serialized version
#    for target_node_name in final_in_nodes:
#         build_IN_levels(target_node_name, dendogram_com, dna_int_graph, genome_fa, out_dir, run_name, json_draft_name, open_peaks_file, picked_tissue, ctcf2bounds, rough_3D_loc)

    #Annotate Global graph, load from file and annotate  
    json_comm = open(json_draft_name,"r")
    json_in = json.load(json_comm)
    json_comm.close()
    draft_genome_out_annot = json_graph.node_link_graph(json_in)
    draft_genome_out_annot.graph["complete"] = False

    #add colored nodes by chr
    with open("all_INs.bed","w") as all_INs:
        for node in draft_genome_out_annot.nodes():
            arr = re.split(r"[-:]",node)
            draft_genome_out_annot.nodes[node]["name"] = str(node) + " " + dna_int_graph.nodes[node]["in_name"]
            all_INs.write(arr[0] + "\t" + str(arr[1]) + "\t" + str(arr[2]) + "\n")

    os.system("sort -u all_INs.bed | sort-bed - > all_INs_sorted.bed")
    os.remove("all_INs.bed")
    color_chrom = {"chr1": "rgb(240,163,255)",
                     "chr2": "rgb(0,117,220)",
                     "chr3": "rgb(153,63,0)",
                     "chr4": "rgb(76,0,92)",
                     "chr5": "rgb(25,25,25)",
                     "chr6": "rgb(0,92,49)",
                     "chr7": "rgb(43,206,72)",
                     "chr8": "rgb(255,204,153)",
                     "chr9": "rgb(128,128,128)",
                     "chr10": "rgb(148,255,181)",
                     "chr11": "rgb(143,124,0)",
                     "chr12": "rgb(157,204,0)",
                     "chr13": "rgb(194,0,136)",
                     "chr14": "rgb(0,51,128)",
                     "chr15": "rgb(255,164,5)",
                     "chr16": "rgb(255,168,187)",
                     "chr17": "rgb(66,102,0)",
                     "chr18": "rgb(255,0,16)",
                     "chr19": "rgb(94,241,242)",
                     "chr20": "rgb(0,153,143)",
                     "chr21": "rgb(224,255,102)",
                     "chr22": "rgb(116,10,255)",
                     "chrX": "rgb(153,0,0)",
                     "chrY": "rgb(255,255,128)",
                     "chrM": "rgb(255,255,0)"}

    sorted_ins_list = list()

    #Write IN area out bed file
    IN_minmax_out = (out_dir + run_name)[0:80]+ "_global_IN.bed"
    with open(IN_minmax_out, "w") as all_INs_bed:
        with open("all_INs_sorted.bed","r") as all_INs_sorted:
            for IN_node in all_INs_sorted:
                #set size of nodes according to sum of open peaks
                node_arr = IN_node.strip().split()
                node_name = (str(node_arr[0]) + ":" + str(node_arr[1]) + "-" + str(node_arr[2]))
                with open(draft_genome_out_annot.nodes[node_name]["json_name"]) as json_IN_file:
                    json_in = json.load(json_IN_file)
                    cur_IN_graph = json_graph.node_link_graph(json_in)
                #write min/max blocks of INs to bed file
                out_in_bed = (str(cur_IN_graph.graph["community_chr"]) + "\t" + 
                              str(cur_IN_graph.graph["community_min"]) + "\t" + 
                              str(cur_IN_graph.graph["community_max"]) + "\t" + 
                              cur_IN_graph.graph["target_gene"] + "\n")
                all_INs_bed.write(out_in_bed)
                #sum size of open regions
                sum_open = 0
                for node in cur_IN_graph.nodes(data=True):
                    if "open_size" in node[1]:
                        sum_open += node[1]["open_size"]
                #Set global graph attributes
                draft_genome_out_annot.nodes[node_name]["open_size"] = sum_open
                draft_genome_out_annot.nodes[node_name]["size"] = 20
                draft_genome_out_annot.nodes[node_name]["color"] = color_chrom[node_arr[0]]
                sorted_ins_list.append(node_name)
            
    #Add chromosome edges connecting INs 
    for i in range(0,len(sorted_ins_list)-2):
        IN1_arr = re.split(r"[-:]",sorted_ins_list[i])
        IN2_arr = re.split(r"[-:]",sorted_ins_list[i+1])
        if IN1_arr[0] == IN2_arr[0]:
            draft_genome_out_annot.add_edge(sorted_ins_list[i],sorted_ins_list[i+1])

    #Output annotated global graph as JSON output
    #Flag as complete
    draft_genome_out_annot.graph["complete"] = True
    json_out_graph = json_graph.node_link_data(draft_genome_out_annot)
    json_out = json.dumps(json_out_graph)
    json_draft_annot_file = run_name[0:80]+ "_global_graph_annotated.json"
    json_draft_annot_name = (out_dir + json_draft_annot_file)
    json_draft_annot = open(json_draft_annot_name,"w+")
    with open(json_draft_annot_name,"w+") as json_draft_annot:
        json_draft_annot.write(json_out)

    #Render global graph to plotly html
    in_json_out_plotly_fig(json_draft_annot_name, 
            make_3D=True, global_in=True, html_out="annotated_global_3D_graph.html")
    
    #Generate 2D genome views:
    in_json_out_genome_views(json_draft_annot_file,
                             out_dir,
                             home_dir,
                             out_dir,
                             run_name
                            )
    
    os.remove("all_nodes_sorted.bed")
    os.remove("open_peaks_nodes.bed")
    os.remove("ctcf_loops.bed")
    os.remove("all_INs_sorted.bed")
