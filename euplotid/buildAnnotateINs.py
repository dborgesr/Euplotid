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
#depends on bedops and deepbind too!

#Function structure:
#Helper / resource functions
#SNP/CNV functions
#Open regions functions
#Insulated Neighborhood functions
#Whole genome functions

def get_target_nodes(G,args,cell_type):
    target_nodes = list()
    if args["top_gene_fpkm"] is not None:
        target_nodes = [node[0] for node in G.nodes(data=True) if "fpkm" in node[1] and node[1]["fpkm"] > float(args["top_gene_fpkm"])]
    elif args["marker_region"] is not None:
        marker_arr = re.split(r"[-:]",args["marker_region"])
        chr_marker = marker_arr[0]
        marker_start = int(marker_arr[1])
        marker_end = int(marker_arr[2])
        for node in G.nodes(data=True):
            arr = re.split(r"[-:]",node[0])
            node_chr = arr[0]
            node_start = int(arr[1])
            node_end = int(arr[2])
            if "fpkm" in node[1] and node_chr == chr_marker and node_start >= marker_start and node_end <= marker_end:
                target_nodes.append(node[0])
    elif args["insulated_neigh"] is not None:
        target_names = args["insulated_neigh"].strip(",").split(",")
        for node in G.nodes(data=True):
            if "name" in node[1]:
                node_names = set(node[1]["name"].split(","))
                if len(node_names.intersection(target_names)) > 0:
                    target_nodes.append(node[0])
    else:
        if cell_type == "primed":
            source_gene_names = "SMAD3,HMGB3,TBX3,LEFTY1,KLF4,NANOG,PRDM14,SOX2,POU5F1,HIF1A,DUSP6,SOX17,IRX3,PAX3"
        elif cell_type == "naive":
            source_gene_names = "SMAD3,HMGB3,TBX3,LEFTY1,KLF4,NANOG,PRDM14,SOX2,POU5F1,HIF1A,DUSP6,SOX17,IRX3,PAX"
        elif cell_type == "neuron":
            source_gene_names = "SNCA,MLH1,PARK2,LRRK2,PINK1"
        elif cell_type == "npc":
            source_gene_names = "SNCA,MLH1,PARK2,LRRK2,PINK1"
        elif cell_type == "gm12878":
            source_gene_names = "LMO2,TAL1,NOTCH1,FGFR1,RUNX1,GATA3"
        elif cell_type == "jurkatt":
            source_gene_names = "LMO2,TAL1,NOTCH1,FGFR1,RUNX1,GATA3"
        elif cell_type == "k562":
            source_gene_names = "TP53,RB1,WT1,NF1,NF2,APC,TSC1,TSC2,DCC,BCRA1,BCRA2,STK11,MSH2,MLH1,VHL,CDKN2A,PTCH1,MEN1,MET,PTEN,ATM,BLM,XPA,ERCC3,XPC,ERCC2,ERCC4,ERCC5,POLH,FANCA,MYC,EGFR,KRAS,HRAS"
        elif cell_type == "hela":
            source_gene_names = "TP53,RB1,WT1,NF1,NF2,APC,TSC1,TSC2,DCC,BCRA1,BCRA2,STK11,MSH2,MLH1,VHL,CDKN2A,PTCH1,MEN1,MET,PTEN,ATM,BLM,XPA,ERCC3,XPC,ERCC2,ERCC4,ERCC5,POLH,FANCA,MYC,EGFR,KRAS,HRAS"
        elif cell_type == "mesc":
            source_gene_names = "nanog,sox2,pou5f1,prdm14,nlrp12,myc,kras,egfr,hras,actb,gapdh"
        else:
            raise NameError("Cant find that cell type, try: primed, naive, neuron, NPC, GM12878, Jurkatt, K562, HeLa, mESC")
        default_gene_names = [x.upper() for x in source_gene_names.split(",")]
        for node in G.nodes(data=True):
            if "name" in node[1]:
                node_names = node[1]["name"].split(",")
                if len(set(node_names).intersection(default_gene_names)) > 0:
                    target_nodes.append(node[0])
    return target_nodes

def liftover_chain(gen_loc, chain_file):
    arr = re.split(r"[-:]",gen_loc)
    lo = LiftOver(chain_file)
    hg38_start = lo.convert_coordinate(arr[0], int(arr[1])-1)[0][1]
    hg38_end = lo.convert_coordinate(arr[0], int(arr[2])-1)[0][1]
    gen_out = arr[0] + ":" + str(hg38_start) + "-" + str(hg38_end)
    return gen_out

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
        
#overlap bed file w/ graph and return list of overlapping nodes
def overlap_graph_bed(G, bed):
    os.system("sort-bed " + bed + " | bedops -e -1 " + "all_nodes_sorted.bed - > " + "overlapped_nodes.bed")
    overlapped_nodes = set()
    with open("overlapped_nodes.bed","r") as over_bed:
        for node in over_bed:
            arr = node.strip().split()
            overlapped_nodes.add(str(arr[0]) + ":" + str(arr[1]) + "-" + str(arr[2]))
    os.remove("overlapped_nodes.bed")
    return overlapped_nodes

def label_gene_nodes(G, bed):
    tss_slop = 50000
    source_in_file = open("source_tss.bed","w")
    gene_names = list()
    fpkm_genes = list()
    #extract TSS from bed file
    for gene in open(bed):
        arr = gene.strip().split()
        if arr[3] == "+" and (arr[0] in chr_names):
            source_in_file.write(arr[0] + "\t" + str(int(arr[1])) + "\t" + str(int(arr[1])+1) + 
                                 "\t" + str(arr[3]) + "\t" + str(arr[4]) + "\t" + str(arr[5]) + "\n")
        elif arr[3] == "-" and (arr[0] in chr_names):
            source_in_file.write(arr[0] + "\t" + str(int(arr[2])) + "\t" + str(int(arr[2])+1) + 
                                 "\t" + str(arr[3]) + "\t" + str(arr[4]) + "\t" + str(arr[5]) + "\n")
        gene_names.append(arr[4].upper())
        fpkm_genes.append(float(arr[5]))
    source_in_file.close()
    os.system("sort-bed source_tss.bed | closest-features --dist --closest - all_nodes_sorted.bed > closest_nodes.bed")
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
                if "name" not in G.node[node_name]:
                    G.node[node_name]["name"] = set()
                G.node[node_name]["gene"] = True
                G.node[node_name]["name"].add(gene_tss_arr[4].upper())
                G.node[node_name]["fpkm"] = float(gene_tss_arr[5])
    #join gene names
    for node in G.nodes(data=True):
        if ("name" in node[1]):
            G.node[node[0]]["name"] = ",".join(G.node[node[0]]["name"])
    #cleanup
    os.remove("source_tss.bed")
    os.remove("closest_nodes.bed")
    return G

def add_TFBS_graph(G, TFBS_bed):
    os.system("bedmap --echo --echo-map-id-uniq --multidelim , --bp-ovr 1 all_nodes_sorted.bed " + TFBS_bed + " > nodes_tfbs.bed")
    with open("nodes_tfbs.bed") as nodes_all_tfbs:
        for node_tf in nodes_all_tfbs:
            arr = node_tf.strip().split("|")
            node_arr = arr[0].split()
            node_name = node_arr[0] + ":" + node_arr[1] + "-" + node_arr[2]
            if len(arr)>1:
                G.node[node_name]["tfbs"] = arr[1].strip()
            else:
                G.node[node_name]["tfbs"] = "NA"
    #cleanup
    os.remove("nodes_tfbs.bed")
    return G

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
            G.node[id_node]["x"] = global_pos[0]
            G.node[id_node]["y"] = global_pos[1]
            G.node[id_node]["z"] = global_pos[2]
    os.remove(sub_nodes)
    os.remove(xyz_nodes)
    return G

def add_variants_predict(G, homo_gen, chain_file, TF_RBP_ids, tissue_type, target_node_name,out_dir):
    #Function to add variants and predict their effect using DeepBind
    in_genes = G.node[target_node_name]["in_name"].split(",")
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
    basset_name_out = os.getcwd() + "/" + target_node_name.replace(":","_").replace("-","_") + "_basset"
    if not os.path.exists(basset_name_out):
        os.makedirs(basset_name_out)
    open_vcf = open(basset_name_out + "/open_variants.vcf","w+")
    #find all variants that fall within these open peaks using MyVariant.info 
    mv = myvariant.MyVariantInfo()
    for node in G.nodes(data=True):
        if "dist2anchor" in node[1]:
            #Pull all variants within chromatin accessibility peak
            all_snps = mv.query(node[0], fetch_all=True)
            for snp in all_snps:
                if "dbsnp" not in snp:
                    continue
                add_any_SNPs = True      
                #Add CNVs/SNVs to graph
                G.add_node(snp["_id"],gen_start=snp["hg19"]["start"], all_ref = str(snp["dbsnp"]["ref"]),
                           all_alt = str(snp["dbsnp"]["alt"]),gen_end=snp["hg19"]["end"], gen_chrom="chr"+str(snp["chrom"]), 
                           name=snp["_id"])
                G.add_edge(snp["_id"], node[0],label="",weight=1.0)
                #print to VCF
                open_vcf.write("chr"+str(snp["chrom"]) + "\t" + str(snp["hg19"]["start"]) + "\t" + snp["_id"] +
                               "\t" + str(snp["dbsnp"]["ref"]) + "\t" + str(snp["dbsnp"]["alt"]) + "\n")
                
                snp_loc = str(G.node[snp["_id"]]["gen_chrom"]) + ":" + str(G.node[snp["_id"]]["gen_start"]) + "-" + str(G.node[snp["_id"]]["gen_end"])
                #add Differentially Methylated Cytosine (DMC) data from NGSMethDB
                if G.graph["species"] == "Human2":
                    DMC_var = get_NGSmethDMC(snp_loc,chain_file,chain_file2)
                    G.node[snp["_id"]]["DMC_node"] = DMC_var
                if ("cadd" in snp):
                    G.node[snp["_id"]]["cadd"] = snp["cadd"]["rawscore"]
                #add rsid if present and add GTEX eQTL pval and if target gene has Ensembl ID
                if (ensembl_target is not None) and ("dbsnp" in snp) and ("rsid" in snp["dbsnp"]):
                    serv = "http://rest.ensembl.org/eqtl/variant_name/homo_sapiens/"
                    ext = snp["dbsnp"]["rsid"].strip() + "?content-type=application/json;statistic=p-value;stable_id=" + ensembl_target + ";tissue=" + tissue_type 
                    r = requests.get(serv + ext, headers={ "Content-Type" : "application/json"})
                    if r.ok:
                        ret_val = r.json()
                        if isinstance(ret_val,list) and (len(ret_val)>0):
                            eqtl_pval = ret_val[0]["value"]
                            G.node[snp["_id"]]["gtex_eqtl_pval"] = round(eqtl_pval,6)
                    G.node[snp["_id"]]["rsid"] = snp["dbsnp"]["rsid"]
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
                    G.node[snp["_id"]]["grasp_pheno"] = pheno_gwas.replace(" ","_")
                    G.node[snp["_id"]]["grasp_pmid"] = pmid_gwas.replace(" ","_")
    open_vcf.close()
    #predict effect of all variation added to graph using DeepBind and Basset
    if add_any_SNPs:
        print("Deepbind SNP prediction of: " + target_node_name)
        G = deepbind_predict_SNPs(G, homo_gen, TF_RBP_ids,target_node_name)
        print("Basset SNP prediction of: " + target_node_name)
        G = add_basset_sad_sat(G, homo_gen, target_node_name,out_dir)
    try:
        os.remove(basset_name_out + "open_variants.vcf")
    except OSError:
        pass
    shutil.rmtree(basset_name_out)
    return G

def deepbind_predict_SNPs(G, homo_gen, TF_RBP_ids,target_node_name):                        
    #predict effect of CNVs/SNPs which were added to Insulated Neighborhood
    flank_snp = 15
    mut_seqs = list()
    ref_seqs = list()
    snp_names = list()
    mut_snps = open(target_node_name + "_mut_snps.fa", "w+")
    ref_snps = open(target_node_name + "_ref_snps.fa", "w+")
    #get deepbind id to TF name mapping
    db_id2name = dict()
    with open(TF_RBP_ids) as db_dp:
        for db_id in db_dp:
            db_id2name[db_id.split("#")[0].strip()] = db_id.split("#")[1].strip()
    for node in G.nodes(data=True):
        if "all_ref" in node[1]:
            chrom = node[1]["gen_chrom"]
            snp_loc = node[1]["gen_start"]
            snp_names.append(node[0])
            #get 5 tiled sequences and average the probability across them
            for offset in range(-2,3):
                seq_start = int((snp_loc-flank_snp)+(offset*5.0))
                seq_end = int((snp_loc+flank_snp+1)+(offset*5.0))
                snp_seq = Seq.Seq(homo_gen.fetch(chrom,seq_start,seq_end))
                ref_seq = snp_seq.tomutable()
                ref_req = SeqRecord(ref_seq, node[0],"","")
                ref_seqs.append(ref_req)
                mut_seq = snp_seq.tomutable()
                alt_seq = node[1]["all_alt"]
                if len(alt_seq) > 1:
                    max_ix = min(len(mut_seq),(flank_snp-1+len(alt_seq)))
                    mut_seq[flank_snp-1:max_ix] = alt_seq
                else:
                    mut_seq[flank_snp-1] = node[1]["all_alt"]
                mut_req = SeqRecord(mut_seq, node[0],"","")
                mut_seqs.append(mut_req)
    #write sorrounding sequences of SNP to file and reopen again for deepbind
    SeqIO.write(mut_seqs, mut_snps, "fasta")
    SeqIO.write(ref_seqs, ref_snps, "fasta")
    mut_snps.close()
    ref_snps.close()
    #call deepbind and dump output
    os.system("/root/deepbind/deepbind " + TF_RBP_ids + " " +  target_node_name + "_mut_snps.fa" + " > " + target_node_name + "_deep_mut.txt")
    os.system("/root/deepbind/deepbind " + TF_RBP_ids + " " +  target_node_name + "_ref_snps.fa" + " > " + target_node_name + "_deep_ref.txt")
    #open deepbind output for ref and mut
    deep_mut_out = open(target_node_name + "_deep_mut.txt","r+")
    deep_ref_out = open(target_node_name + "_deep_ref.txt","r+")
    #peel header off
    header = deep_mut_out.readline().strip().split() 
    header = deep_ref_out.readline().strip().split() 
    #suck deebind output into numpy array
    deep_mut_prob = np.loadtxt(deep_mut_out,  ndmin = 2)
    deep_ref_prob = np.loadtxt(deep_ref_out,  ndmin = 2)
    deep_mut_out.close()
    deep_ref_out.close()
    for snp_num in range(0,len(snp_names)):
        #compare scores of TF binding according to Deepbind paper
        deep_ref_prob_cut = np.mean(deep_ref_prob[(snp_num*5):((snp_num+1)*5),:],axis=0)
        deep_mut_prob_cut = np.mean(deep_mut_prob[(snp_num*5):((snp_num+1)*5),:],axis=0)
        max_diffs = np.maximum.reduce([np.zeros(len(deep_ref_prob_cut)),deep_ref_prob_cut,deep_mut_prob_cut])      
        deep_cut_prob = (deep_ref_prob_cut - deep_mut_prob_cut) * max_diffs
        max_shown = min(5, sum(abs(deep_cut_prob)>=3))
        if max_shown > 0:
            temp_sort = np.array(abs(deep_cut_prob)).argsort()[::-1]
            deep_probs_ix_sort = [ix for ix in temp_sort if np.isfinite(deep_cut_prob[ix])][0:max_shown]
            out_prob_name = ""
            top_tf_names = [db_id2name[header[ix]] for ix in deep_probs_ix_sort]
            top_tf_probs = [deep_cut_prob[ix] for ix in deep_probs_ix_sort]
            for x in zip(top_tf_names,top_tf_probs):
                out_prob_name += str(x[0]) + ":" + str(round(x[1],2)) + ", "
            G.node[snp_names[snp_num]]["deep_score"] = out_prob_name
    os.remove(target_node_name + "_deep_mut.txt")
    os.remove(target_node_name + "_deep_ref.txt")
    os.remove(target_node_name + "_mut_snps.fa")
    os.remove(target_node_name + "_ref_snps.fa")
    return G

def add_basset_sad_sat(G, homo_gen, target_node_name,out_dir):
    basset_out_name = out_dir + "/" + target_node_name.replace(":","_").replace("-","_") + "_basset"
    target_in_name = dna_int_graph.node[target_node_name]["in_name"]
    #trained DNAse peak model
    model_file = "/input_dir/pretrained_model.th"
    #Sequences stored in HDF5 format of encode DNAse peaks
    seqs_file = "/input_dir/encode_roadmap.h5"
    #table of DNAse target BED files used to train model
    targets_file = "/root/Basset/tutorials/sad_eg/sample_beds.txt"
    targetID_celltype = dict()
    count=0
    with open(targets_file) as id2type:
        for c_type in id2type:
            arr = c_type.strip().split()
            targetID_celltype[arr[0]] = count
            count += 1
    cmd = ("/root/Basset/src/basset_sad.py -f " + (homo_gen.filename).decode("utf-8") + " -l 600 -o " + basset_out_name + " -t " + targets_file + " " +  model_file + " " + basset_out_name + "/open_variants.vcf")
    print(cmd)
    os.system(cmd)
    #Pick SNP maximizing for sum of Delta SAD across all cell types in trained model
    sad_table = pd.read_table(basset_out_name + "/sad_table.txt", delim_whitespace=True, header = 0)
    sad_table_dense = sad_table.pivot_table(index='rsid', columns='target', values='pred')
    serial_sad_table = pickle.dumps(sad_table_dense, protocol=0)
    G.graph["sad_table"] = serial_sad_table
    sad_table["sad_abs"] = abs(sad_table["pred"])
    sad_rsid_table = sad_table.groupby("rsid").sum()
    sad_rsid_table = sad_rsid_table.sort_values(by="sad_abs",ascending=False)
    #Store SAD in each SNP node
    top_snp_names = sad_rsid_table.index
    for snp_index, snp_row in sad_rsid_table.iterrows():
        G.node[snp_index]["sad_abs_sum"] = snp_row["sad_abs"]
        G.node[snp_index]["sad_pred"] = ""
    for sad_index, sad_row in sad_table.iterrows():
        G.node[sad_row["rsid"]]["sad_pred"] += str(round(sad_row["pred"],4)) + ","
    target_sad_name = (top_snp_names[0].replace(":","_")).replace(">","_")
    top_snp_sad_abs_sum = sad_rsid_table["sad_abs"] 
    top_snp_sad_profile = sad_table[sad_table["rsid"].str.contains(top_snp_names[0])]
    top_snp_ix = top_snp_sad_profile["sad_abs"].argmax()
    top_snp_sad_table = top_snp_sad_profile.ix[top_snp_ix] 
    top_snp_ctype_target = targetID_celltype[top_snp_sad_table["target"]]
    top_snp_node = G.node[top_snp_names[0]]
    with open(basset_out_name + "/top_variant.vcf","w") as filt_out:
        filt_out.write(str(top_snp_node["gen_chrom"]) + "\t" + str(top_snp_node["gen_start"]) + "\t" + top_snp_names[0] + 
                        "\t" + top_snp_node["all_ref"] + "\t" + top_snp_node["all_alt"] + "\n")
    #flag top SNP in open peak
    G.node[top_snp_names[0]]["top_open_snp"] = True
    #Perform in-silico mutagenesis on top SNP to generate output PDFs
    cmd = ("/root/Basset/src/basset_sat_vcf.py -f " + (homo_gen.filename).decode("utf-8") + " -t " + str(top_snp_ctype_target) + " -o " + basset_out_name + 
           " " + model_file + " " + basset_out_name + "/top_variant.vcf")
    os.system(cmd)
    #Save output png & cleanup
    #tag in JSON file too
    G.node[top_snp_names[0]]["sad_pdf"] = target_in_name + "_sad_heat.png"
    out_pdfs = [pdf for pdf in glob.glob(basset_out_name+"/*heat.png") if target_sad_name in pdf]
    label_snp = ["ref","alt"]
    count = 1
    for pdf in out_pdfs:
        pdf_sat = pdf.split("/")
        G.node[top_snp_names[0]]["sad_mut_"+label_snp[count%2]] = str(pdf_sat[len(pdf_sat)-1])
        count += 1
        if os.path.isfile(pdf):
            os.rename(pdf, out_dir+"/"+pdf_sat[len(pdf_sat)-1])
    return G

def add_openRegions_predict(G, homo_gen, chain_file, chain_file2, TF_RBP_ids,open_peaks, target_node_name, in_level, openpeak2bassetpeak):
    #Function to add open chromatin accessibility peaks to graph and predict TF binding using DeepBind
    #how far can the chromatin accessibility be from the node to be considered?
    open_region_slop = 100000
    orig_edges = G.edges()
    with open(open_peaks,"r") as all_peaks:
        for peak in all_peaks:
            arr = peak.strip().split("|")
            if len(arr) > 1 and (arr[1] != "NA") and (abs(int(arr[2])) <= open_region_slop):
                dist2anchor = abs(int(arr[2]))
                arr_atac = arr[0].strip().split()
                node_arr = arr[1].strip().split()
                atac_node = arr_atac[0] + ":" + str(int(arr_atac[1])) + "-" + str(int(arr_atac[2]))
                open_size = float(arr_atac[3])
                anchor_node = node_arr[0] + ":" + str(int(node_arr[1])) + "-" + str(int(node_arr[2]))
                if anchor_node in G.nodes():
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
                    #add Differentially Methylated Cytosine (DMC) data from NGSMethDB if human
                    if G.graph["species"] == "Human2":
                        DMC_node = get_NGSmethDMC(atac_node,chain_file,chain_file2)
                    else:
                        DMC_node = ""
                    top_tf_names=("|".join(top_tf_names))
                    top_tf_probs=("|".join(top_tf_probs))
                    filt_tf_open = top_tf_names[0:min([len(top_tf_names),15])]
                    G.add_node(atac_node, mid=atac_mid,name=atac_node, dist2anchor=dist2anchor,
                              deepbind_tf=predict_tf_deepbind,top_tf_names=top_tf_names, filt_tf_open = filt_tf_open,
                               top_tf_probs=top_tf_probs, DMC_node=DMC_node, open_size=open_size)
                    G.add_edge(anchor_node,atac_node,label="",weight=1.0)
    os.system("bedmap --echo --echo-map-id-uniq --multidelim , --bp-ovr 1 " + open_peaks + " " + encode_bed + " > " + target_node_name + "_tfbs.bed")
    with open(target_node_name + "_tfbs.bed","r") as open_tfbs_all:
        for node_tf in open_tfbs_all:
            arr = node_tf.strip().split("|")
            node_arr = arr[0].split()
            node_name = node_arr[0] + ":" + node_arr[1] + "-" + node_arr[2]
            if len(arr)>1 and node_name in G:
                G.node[node_name]["tfbs"] = arr[1].strip()
    #cleanup
    os.remove(target_node_name + "_tfbs.bed")
    
    #add predicted co-occurences
    for edge in orig_edges:
        left_tf = set()
        right_tf = set()
        left_tf_enc = set()
        right_tf_enc = set()
        #overlap top 15 TFs of each chromatin accessibility site as ordered by deepbind score
        for neigh in G[edge[0]]:
            if "top_tf_names" in G.node[neigh]:
                left_tf |= set(G.node[neigh]["top_tf_names"].split("|")[0:15])
            if "tfbs" in G.node[neigh] and "top_tf_names" in G.node[neigh]:
                left_tf_enc |= set(G.node[edge[0]]["tfbs"].strip().split(","))
        for neigh in G[edge[1]]:
            if "top_tf_names" in G.node[neigh]:
                left_tf |= set(G.node[neigh]["top_tf_names"].split("|")[0:15])
            if "tfbs" in G.node[neigh] and "top_tf_names" in G.node[neigh]:
                right_tf_enc |= set(G.node[edge[0]]["tfbs"].strip().split(","))
        overlapped_tf = left_tf.intersection(right_tf)
        overlapped_tf_enc = left_tf_enc.intersection(right_tf_enc)
        if len(overlapped_tf) > 0:
            over_names = ""
            for tf in overlapped_tf:
                over_names += tf + ", "
            G[edge[0]][edge[1]]["overlapped_tf"] = over_names
        
        if len(overlapped_tf_enc) > 0:
            over_names_enc = ""
            for tf in overlapped_tf_enc:
                over_names_enc += tf + ", "
            G[edge[0]][edge[1]]["overlapped_tf_enc"] = over_names_enc
        
        if len(left_tf)>0:
            G.node[edge[0]]["filt_tf_peak"] = "|".join(left_tf)
        if len(right_tf)>0:
            G.node[edge[1]]["filt_tf_peak"] = "|".join(right_tf)
            
        if len(left_tf_enc)>0:
            G.node[edge[0]]["tfbs_enc_all"] = "|".join(left_tf_enc)
        if len(right_tf_enc)>0:
            G.node[edge[1]]["tfbs_enc_all"] = "|".join(right_tf_enc)
            
    #remove long prediction of annotations
    for node in G.nodes():
        if "top_tf_names" in G.node[node]:
            G.node[node]["top_tf_names"] = ""
            G.node[node]["top_tf_probs"] = ""
    return G

def deepbind_predict_tf_range(chrom, start, end, homo_gen, TF_RBP_ids, target_node_name):
    size_window = 30
    stagger_window = 5
    deep_node_name = target_node_name.replace(":","_").replace("-","_")
    #filt TFs out w/ scores below this threshold
    tf_filt = 2.0
    open_seqs = list()
    for left_bound in range(int(start),int(end),stagger_window):
        open_seq = Seq.Seq(homo_gen.fetch(chrom,left_bound,left_bound+size_window))
        ref_req = SeqRecord(open_seq, chrom+":"+str(left_bound)+"-"+str(left_bound+size_window),"","")
        open_seqs.append(ref_req)
    open_seqs_fa = open(deep_node_name + "_open_anchor.fa", "w+")
    #get deepbind id to TF name mapping
    db_id2name = dict()
    with open(TF_RBP_ids) as db_dp:
        for db_id in db_dp:
            db_id2name[db_id.split("#")[0].strip()] = db_id.split("#")[1].strip()
    SeqIO.write(open_seqs, open_seqs_fa, "fasta")
    open_seqs_fa.close()
    #call deepbind and dump output
    os.system("/root/deepbind/deepbind " + TF_RBP_ids + " " +  deep_node_name + "_open_anchor.fa" + " > " + deep_node_name + "_deep_open.txt")
    deep_open_out = open(deep_node_name + "_deep_open.txt","r+")
    header = deep_open_out.readline().strip().split()
    deep_open_prob = np.loadtxt(deep_open_out,  ndmin = 2)
    deep_open_out.close()
    #smooth over overlapping windows and flatten array
    deep_open_prob_smooth = list()
    for tf in range(0,np.shape(deep_open_prob)[1]):
        #deep_open_prob_smooth.extend(moving_average(deep_open_prob[:,tf],n=4))
        deep_open_prob_smooth.extend(pd.rolling_max(deep_open_prob[:,tf],4)[3:len(deep_open_prob[:,tf])])
    deep_probs_ix_sort = np.array(deep_open_prob_smooth).argsort()[::-1]
    tf_names = [db_id2name[header[ix%np.shape(deep_open_prob)[1]]] 
                for ix in deep_probs_ix_sort if deep_open_prob_smooth[ix] >= tf_filt]
    tf_probs = [deep_open_prob_smooth[ix] for ix in deep_probs_ix_sort 
                if deep_open_prob_smooth[ix] >= tf_filt]    
    tf_loc =  [(np.floor(ix/np.shape(deep_open_prob)[1]) + int(start)) 
               for ix in deep_probs_ix_sort if deep_open_prob_smooth[ix] >= tf_filt]
    os.remove(deep_node_name + "_deep_open.txt")
    os.remove(deep_node_name + "_open_anchor.fa")
    return [tf_names,tf_probs,tf_loc]

def basset_predict_tf_range(open_region, homo_gen, target_node_name, openpeak2bassetpeak):
    #trained DNAse peak model
    model_file = "/input_dir/pretrained_model.th"
    #Sequences stored in HDF5 format of encode DNAse peaks
    seqs_file = "/input_dir/encode_roadmap.h5"
    #table of DNAse target BED files used to train model
    targets_file = "/root/Basset/tutorials/sad_eg/sample_beds.txt"
    #expand region to at least 600 bp, default size model was trained on
    homo_gen = pysam.FastaFile("/input_dir/hg19.fa")
    arr = re.split(r"[-:]",open_region)
    chrom = arr[0]
    start = arr[1]
    end = arr[2]
    stagger_window = 10
    size_window = 600
    open_seqs = list()
    for left_bound in range(int(start),int(end),stagger_window):
        open_seq = Seq.Seq(homo_gen.fetch(chrom,left_bound,left_bound+size_window))
        ref_req = SeqRecord(open_seq, chrom+":"+str(left_bound)+"-"+str(left_bound+size_window),"","")
        open_seqs.append(ref_req)
    with open("basset_open_anchor.fa", "w+") as open_seqs_fa:
        SeqIO.write(open_seqs, open_seqs_fa, "fasta")
    cmd = ("/root/Basset/src/seq_hdf5.py -r -c -v " + str(len(open_seqs)) + " -t " + str(len(open_seqs)) 
           + " basset_open_anchor.fa /input_dir/encode_roadmap_act.txt open_region.h5")
    subprocess.call(cmd, shell=True)
    cmd = ("/root/Basset/src/basset_motifs.py -s " + str(len(open_seqs)) + " -t -o motifs_out + " 
           + model_file + " open_region.h5")
    subprocess.call(cmd, shell=True)

def get_NGSmethDMC(open_region,chain_file,chain_file2):
    open_region_hg38 =  liftover_chain(open_region,chain_file)
    if open_region_hg38 is None:
        return ""
    ngs_methdb = "http://bioinfo2.ugr.es:8888/NGSmethAPI/hg38/" + open_region_hg38
    meth_db_http = requests.get(ngs_methdb)
    meth_db_out = ""
    if meth_db_http.ok:
        out_meth = meth_db_http.json()
        if len(out_meth)>0:
            for loc in out_meth:
                if "diffmeth_cg" in loc:
                    tissues_meth = set()
                    pval_meth = list()
                    for samp_a in loc["diffmeth_cg"]:
                        for samp_b in loc["diffmeth_cg"][samp_a]:
                            if "methylKit" in loc["diffmeth_cg"][samp_a][samp_b] and (float(loc["diffmeth_cg"][samp_a][samp_b]["methylKit"])<.01):
                                pval_meth.append(float(loc["diffmeth_cg"][samp_a][samp_b]["methylKit"]))
                                tissues_meth |= set(samp_a.split("#"))
                                tissues_meth |= set(samp_b.split("#"))
                            if "MOABS_sim" in loc["diffmeth_cg"][samp_a][samp_b] and (float(loc["diffmeth_cg"][samp_a][samp_b]["MOABS_sim"])<.01):
                                pval_meth.append(float(loc["diffmeth_cg"][samp_a][samp_b]["MOABS_sim"]))
                                tissues_meth |= set(samp_a.split("#"))
                                tissues_meth |= set(samp_b.split("#"))
                    if len(tissues_meth)>0:
                        arr = re.split(r"[-:]",open_region)
                        hg38_pos = arr[0] + ":" + str(loc["pos"]) + "-" + str(loc["pos"]+1)
                        hg19_pos = liftover_chain(hg38_pos,chain_file)
                        meth_db_out += "DMC:(" + hg19_pos + ") number_tissues:" + str(len(tissues_meth)) \
                        + " min_pval:" + str(round(min(pval_meth),9)) + "\n"
    return meth_db_out

def annotate_IN(target_node_name, dendogram_com, level, dna_int_graph, add_open, add_snp, homo_gen, chain_file, chain_file2, TF_RBP_ids,
                open_peaks_file, picked_tissue, ctcf2bounds, comm_open_out, comm_var_out, comm_nodes_out, comm_edges_out, comm_edges_out_bed, ucsc_session_in,out_dir):    
    #fetch Insulated Neighborhood and annotate
    if level <= 1:
        comm_out_path = out_dir + run_name + "_" + dna_int_graph.node[target_node_name]["in_name"].split(":")[0][0:80] + "_lvl"+str(level+1)
    else:
        comm_out_path = out_dir + run_name + "_" + dna_int_graph.node[target_node_name]["in_name"].split(":")[0][0:80] + "_Top"
    cur_neighborhoods = community.partition_at_level(dendogram_com, level)
    cur_nodes = [nodes for nodes in cur_neighborhoods.keys() if cur_neighborhoods[nodes] == cur_neighborhoods[target_node_name]]
    source_community_graph = dna_int_graph.subgraph(cur_nodes)
    target_name = dna_int_graph.node[target_node_name]["in_name"].split(":")[0]
    source_community_graph.graph["target_gene"] = target_name
    
    orig_edges = source_community_graph.edges()
    if add_open:
        print("Annotating Open Regions")
        source_community_graph = add_openRegions_predict(source_community_graph, homo_gen, chain_file, chain_file2, TF_RBP_ids, open_peaks_file, target_node_name, level, dict())
        sys.stdout.flush()
        if add_snp:
            print("Annotating SNPs/CNVs")
            source_community_graph = add_variants_predict(source_community_graph, homo_gen, chain_file, TF_RBP_ids, picked_tissue, target_node_name,out_dir)
            sys.stdout.flush()
    
    #Color Insulated Neighborhood according to node attributes, make mid distance to target node

    for node in source_community_graph.nodes(data=True):
        source_community_graph.node[node[0]]["style"] = "filled"
        if (node[0] in target_nodes):
            source_community_graph.node[node[0]]["color"] = "rgb(17,166,216)"
            source_community_graph.node[node[0]]["size"] = 50
        elif ("gene" in node[1]):
            source_community_graph.node[node[0]]["color"] = "rgb(13,59,224)"
            source_community_graph.node[node[0]]["size"] = 25
        elif node[0] in enh_nodes:
            source_community_graph.node[node[0]]["color"] = "rgb(242,12,50)"
            source_community_graph.node[node[0]]["size"] = 25
        elif node[0] in ctcf_nodes:
            source_community_graph.node[node[0]]["color"] = "rgb(6,188,3)"
            source_community_graph.node[node[0]]["size"] = 25
        elif "top_open_snp" in node[1]:
            source_community_graph.node[node[0]]["color"] = "rgb(255,0,233)"
            source_community_graph.node[node[0]]["size"] = 12
        elif "dist2anchor" in node[1]:
            source_community_graph.node[node[0]]["color"] = "rgb(170,15,226)"
            source_community_graph.node[node[0]]["size"] = 15
        elif "all_ref" in node[1]:
            source_community_graph.node[node[0]]["color"] = "rgb(232,120,9)"
            source_community_graph.node[node[0]]["size"] = 10
        else:
            source_community_graph.node[node[0]]["color"] = "rgb(174,183,180)"
            source_community_graph.node[node[0]]["size"] = 20
            
    #get x,y,z local scatter_graphcoordinates and add to node
    position=nx.spring_layout(source_community_graph,dim=3)    
    for node in position:
        source_community_graph.node[node]["x"] = str(round(position[node][0],5))
        source_community_graph.node[node]["y"] = str(round(position[node][1],5))
        source_community_graph.node[node]["z"] = str(round(position[node][2],5))
    #Add UCSC session
    source_community_graph.graph["ucsc_session"] = ucsc_session_in
    #JSON output
    json_out_graph = json_graph.node_link_data(source_community_graph)
    json_comm = open(comm_out_path + "_graph.json","w+")
    json_out = json.dumps(json_out_graph)
    json_comm.write(json_out)
    json_comm.close()
    
    #2D Genome Browswer output
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

def build_IN_levels(target_node_name, dendogram_com, dna_int_graph, genome_fa, chain_file, chain_file2, out_dir, run_name, json_draft_name, TF_RBP_ids, open_peaks_file,
                    picked_tissue, ctcf2bounds,hsa_3D_loc,ucsc_session):
    #output to genome browser
    comm_edges_out = open(out_dir + run_name + "_target_community_edges_lvl_1.washu","a+")
    comm_edges_out_lvl1 = open(out_dir + run_name + "_target_community_edges_lvl_2.washu","a+")
    comm_edges_out_lvl2 = open(out_dir + run_name + "_target_community_edges_lvl_Top.washu","a+")
    
    comm_edges_out_bed = open(out_dir + run_name + "_target_community_edges_lvl_1.bed","a+")
    comm_edges_out_lvl1_bed = open(out_dir + run_name + "_target_community_edges_lvl_2.bed","a+")
    comm_edges_out_lvl2_bed = open(out_dir + run_name + "_target_community_edges_lvl_Top.bed","a+")
    
    comm_nodes_out = open(out_dir + run_name + "_target_community_nodes_lvl_1.bed","a+")
    comm_nodes_out_lvl1 = open(out_dir + run_name + "_target_community_nodes_lvl_2.bed","a+")
    comm_nodes_out_lvl2 = open(out_dir + run_name + "_target_community_nodes_lvl_Top.bed","a+")
    
    comm_var_out = open(out_dir + run_name + "_target_community_variants_lvl_1.bed","a+")
    comm_var_out_lvl1 = open(out_dir + run_name + "_target_community_variants_lvl_2.bed","a+")
    comm_var_out_lvl2 = open(out_dir + run_name + "_target_community_variants_lvl_Top.bed","a+")

    comm_open_out = open(out_dir + run_name + "_target_community_open_regions_lvl_1.bed","a+")
    comm_open_out_lvl1 = open(out_dir + run_name + "_target_community_open_regions_lvl_2.bed","a+")
    comm_open_out_lvl2 = open(out_dir + run_name + "_target_community_open_regions_lvl_Top.bed","a+")


    target_name = dna_int_graph.node[target_node_name]["in_name"].split(":")[0]
    chr_in = re.split(r"[-:]",target_node_name)[0]
    ucsc_session_in = ucsc_session + "&position=" + chr_in + ":" + str(dna_int_graph.node[target_node_name]["in_max"]) + "-" + str(dna_int_graph.node[target_node_name]["in_min"])
    
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
                                                  homo_gen, chain_file, chain_file2, TF_RBP_ids, open_peaks_file, picked_tissue, ctcf2bounds, comm_open_out,
                                                  comm_var_out, comm_nodes_out, comm_edges_out, comm_edges_out_bed, ucsc_session_in,out_dir)
    #level 1 
    if len(dendogram_com) > 1:
        level = 1
        add_open = False
        add_snp = False
        source_community_graph_lvl1 = annotate_IN(target_node_name, dendogram_com, level, dna_int_graph, add_open, add_snp,
                                                  homo_gen, chain_file, chain_file2, TF_RBP_ids, open_peaks_file, picked_tissue, ctcf2bounds, comm_open_out_lvl1,
                                                  comm_var_out_lvl1, comm_nodes_out_lvl1, comm_edges_out_lvl1, comm_edges_out_lvl1_bed, ucsc_session_in,out_dir)
    #top level
    if len(dendogram_com) > 2:
        level = len(dendogram_com)-1
        add_open = False
        add_snp = False
        source_community_graph_lvl2 = annotate_IN(target_node_name, dendogram_com, level, dna_int_graph, add_open, add_snp,
                                                  homo_gen, chain_file, chain_file2, TF_RBP_ids, open_peaks_file, picked_tissue, ctcf2bounds, comm_open_out_lvl2,
                                                  comm_var_out_lvl2, comm_nodes_out_lvl2, comm_edges_out_lvl2, comm_edges_out_lvl2_bed, ucsc_session_in,out_dir)
    #Load stored draft, add node, and output, dont add chrX,Y,M,22 (xyz coordinates for chr22 funky)
    comm_out_path = run_name + "_" + dna_int_graph.node[target_node_name]["in_name"].split(":")[0][0:80] + "_lvl1"
    if "chrX" not in target_node_name and "chrY" not in target_node_name and "chr22" not in target_node_name and "chrM" not in target_node_name:
        with open(json_draft_name,"r") as json_comm:
            json_in = json.load(json_comm)
        draft_genome_out = json_graph.node_link_graph(json_in)
        draft_genome_out.add_node(target_node_name, name=target_name, json_name=(comm_out_path + "_graph.json"), edges_ucsc_lvl0 = run_name + "_target_community_edges_lvl_1.bed", nodes_lvl0_bed = run_name+"_target_commmunity_nodes_lvl1.bed",
	 open_lvl0_bed = run_name + "_target_community_open_lvl_1.bed", var_lvl0_bed = run_name + "_target_community_variants_lvl_1.bed")
        #add HSA derived XYZ coordinates for IN collections
        draft_genome_out = add_xyz_loc(draft_genome_out,hsa_3D_loc)
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
    

if __name__ == "__main__":
    #Handle command line arguments
    parser = argparse.ArgumentParser(description="Pipeline to output graphical genomic model from 2D data")
    parser.add_argument("-i","--input_directory", help="Where all the input data lives", required=True)
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
    run_name = args["name_output"]
    out_dir = args["output_directory"] + "/" + cell_type + "/" + run_name + "/"
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    #switch to output directory for analysis
    os.chdir(out_dir)

    #Set up environment variables for Basset
    basset_dir = "/root/Basset"
    os.environ["BASSETDIR"] = basset_dir
    os.environ["PATH"] = os.environ["BASSETDIR"] + "/src:" + os.environ["PATH"]
    os.environ["PYTHONPATH"] = os.environ["BASSETDIR"] + "/src:" + ':'.join(sys.path)
    os.environ["LUA_PATH"] = os.environ["BASSETDIR"] + "/src/?.lua;/opt/conda/bin/lua"

    input_dir = dict()
    chain_file2 = home_dir + "hg38ToHg19.over.chain.gz"
    if cell_type == "primed":
        genome_version = "GRCh37.p13"
        species = "Human"
        ucsc_session = "http://genome.ucsc.edu/cgi-bin/hgTracks?hgS_doOtherUser=submit&hgS_otherUserName=euplotidUCSC&hgS_otherUserSessionName=primed"
        chain_file = home_dir + "hg19ToHg38.over.chain.gz"
        hsa_3D_loc = home_dir + "hES_hsa_400kb_xyz.bed"
        dna_ints = home_dir + "primed_.7_origami.bedpe" #Primed-SMC1-petcount-filtered.bedpe,Primed-SMC1-petcount.bedpe, primed_n2.bedpe
        enh_peaks = home_dir + "H1_H3K27ac_primed_MERGED_peaks.bed" #H1_H3K27ac_primed_MERGED_peaks.bed
        prom_peaks = home_dir + "primed_h3k4me3_peaks.bed" #primed_h3k4me3_peaks.bed
        ctcf_peaks = home_dir + "ENCFF001USS_h1_CTCF_hg19.bed" #ENCFF001USS_h1_CTCF_hg19.bed
        open_peaks = home_dir + "GSM2257291_ATAC1+2_IDR_hg19.bed" #GSM816632_H1_dnase_sorted.bed,H1_atac_sorted.bed,GSM2257291_ATAC1+2_IDR_hg19.bed
        rna_seq = home_dir + "ENCFF000DJY_hg19_tss_ucsc_rnaseq_h1.bed" #ENCFF000DJY_hg19_tss_ucsc_rnaseq_h1.bed
        encode_bed = home_dir + "hg19_encode_TFBSv2.bed"
        genome_fa =  home_dir + "hg19.fa"
        TF_RBP_ids = home_dir + "TF_RBP_human.ids"
        add_snp = True
    elif cell_type == "naive":
        genome_version = "GRCh37.p13"
        species = "Human"
        hsa_3D_loc = home_dir + "hES_hsa_400kb_xyz.bed"
        ucsc_session = "http://genome.ucsc.edu/cgi-bin/hgTracks?hgS_doOtherUser=submit&hgS_otherUserName=euplotidUCSC&hgS_otherUserSessionName=primed"
        chain_file = home_dir + "hg19ToHg38.over.chain.gz"
        dna_ints = home_dir + "Naive-SMC1-petcount-filtered.bedpe" #Naive-SMC1-petcount-filtered.bedpe, naive_n2.bedpe
        enh_peaks = home_dir + "H1_H3K27ac_naive_MERGED_peaks.bed" #H1_H3K27ac_naive_MERGED_peaks.bed
        prom_peaks = home_dir + "naive_h3k4me3_peaks.bed" #naive_h3k4me3_peaks.bed
        ctcf_peaks = home_dir + "naive_CTCF_peaks_xiong_paper.bed" #naive_CTCF_peaks_xiong_paper.bed
        open_peaks = home_dir + "GSM2257291_ATAC1+2_IDR_hg19.bed" #GSM816632_H1_dnase_sorted.bed,H1_atac_sorted.bed,GSM2257291_ATAC1+2_IDR_hg19.bed
        rna_seq = home_dir + "naive_rnaseq_RSEM_nonpolya_ucsc.bed" #naive_rnaseq_RSEM_nonpolya_ucsc.bed
        encode_bed = home_dir + "hg19_encode_TFBSv2.bed"
        genome_fa =  home_dir + "hg19.fa"
        TF_RBP_ids = home_dir + "TF_RBP_human.ids"
        add_snp = True
    elif cell_type == "neuron":
        genome_version = "GRCh37.p13"
        species = "Human"
        ucsc_session = "http://genome.ucsc.edu/cgi-bin/hgTracks?hgS_doOtherUser=submit&hgS_otherUserName=euplotidUCSC&hgS_otherUserSessionName=primed"
        chain_file = home_dir + "hg19ToHg38.over.chain.gz"
        hsa_3D_loc = home_dir + "hES_hsa_400kb_xyz.bed"
        dna_ints = home_dir + "Neuron-SMC1-petcount-filtered.bedpe" #Neuron-SMC1-petcount-filtered.bedpe, Neuron_rep2_hg19.bedpe,Neuron_SMC1_n2.bedpe,Neuron-SMC1-petcount-all.bedpe
        enh_peaks = home_dir + "GSM1831746_neuron_h3k27ac.bed" #GSM1831746_neuron_h3k27ac.bed
        prom_peaks = home_dir + "ENCFF075LEC_NPC_h3k4me3.bed" #ENCFF075LEC_NPC_h3k4me3.bed
        ctcf_peaks = home_dir + "bing_NPC_CTCF_2rpm.bed" #bing_NPC_CTCF_2rpm.bed
        open_peaks = home_dir + "GSM2065328_neuronATAC_allreps.bed" #ENCFF001UUI_hg19_cerebllumn_dnase.bed,GSM2065328_neuronATAC_allreps.bed
        rna_seq = home_dir + "neuron_barres_lab_hg19.bed" #neuron_barres_lab_hg19.bed
        encode_bed = home_dir + "hg19_encode_TFBSv2.bed"
        genome_fa =  home_dir + "hg19.fa"
        TF_RBP_ids = home_dir + "TF_RBP_human.ids"
        add_snp = True
    elif cell_type == "npc":
        genome_version = "GRCh37.p13"
        species = "Human"
        hsa_3D_loc = home_dir + "hES_hsa_400kb_xyz.bed"
        ucsc_session = "http://genome.ucsc.edu/cgi-bin/hgTracks?hgS_doOtherUser=submit&hgS_otherUserName=euplotidUCSC&hgS_otherUserSessionName=primed"
        chain_file = "home_dir" + "hg19ToHg38.over.chain.gz"
        dna_ints = home_dir + "NPC-SMC1-petcount-filtered.bedpe" #NPC-SMC1-petcount-filtered.bedpe,NPC_rep2_hg19.bedpe
        enh_peaks = home_dir + "ENCFF704EQR_NPC_h3k27ac.bed" #ENCFF704EQR_NPC_h3k27ac.bed
        prom_peaks = home_dir + "ENCFF075LEC_NPC_h3k4me3.bed" #ENCFF075LEC_NPC_h3k4me3.bed
        ctcf_peaks = home_dir + "bing_NPC_CTCF_2rpm.bed" #bing_NPC_CTCF_2rpm.bed
        open_peaks = home_dir + "GSM2065328_neuronATAC_allreps.bed" #ENCFF001UUI_hg19_cerebllumn_dnase.bed,GSM2065328_neuronATAC_allreps.bed
        rna_seq = home_dir + "neuron_barres_lab_hg19.bed" #ENCFF947FTB_h9_neuronal_RAMPAGE_fpkm.bed,neuron_barres_lab_hg19.bed
        encode_bed = home_dir + "hg19_encode_TFBSv2.bed"
        genome_fa =  home_dir + "hg19.fa"
        TF_RBP_ids = home_dir + "TF_RBP_human.ids"
        add_snp = True
    elif cell_type == "gm12878":
        genome_version = "GRCh37.p13"
        species = "Human"
        ucsc_session = "http://genome.ucsc.edu/cgi-bin/hgTracks?hgS_doOtherUser=submit&hgS_otherUserName=euplotidUCSC&hgS_otherUserSessionName=primed"
        chain_file = home_dir + "hg19ToHg38.over.chain.gz"
        hsa_3D_loc = home_dir + "hES_hsa_400kb_xyz.bed"
        #GM12878_SMC_HiChip.bedpe, GM12878-CTCF-filtered.bedpe,GM12878_merged_.1_filt_mango_hg19.bedpe
        #ENCFF002EMO_GM12878_rep1_Rad21_Chia_hg19.bedpe,GSE63525_GM12878_primary+replicate_HiCCUPS_looplist.bed,
        #GSE63525_KBM7_HiCCUPS_looplist.bed
        dna_ints = home_dir + "GM12878_SMC_HiChip.bedpe" 
        enh_peaks = home_dir + "ENCFF001SUG_GM12878_H3k27ac_hg19_peaks.bed" #ENCFF001SUG_GM12878_H3k27ac_hg19_peaks.bed
        prom_peaks = home_dir + "ENCFF001SUF_GM12878_H3k4me3_hg19_peaks.bed" #ENCFF001SUF_GM12878_H3k4me3_hg19_peaks.bed
        ctcf_peaks = home_dir + "GM12878-CTCF-filtered.bed" #GM12878-CTCF-filtered.bed
        open_peaks = home_dir + "GSE47753_GM12878_ATACseq_50k_AllReps_readcnts.bed" #GSE66386_GM12878_nucpos.bed,GSE47753_GM12878_ATACseq_50k_AllReps_readcnts.bed
        rna_seq = home_dir + "ENCFF000DAR_GM12878_rnaseq_hg19.bed" #ENCFF000DAR_GM12878_rnaseq_hg19.bed
        encode_bed = home_dir + "hg19_encode_TFBSv2.bed"
        genome_fa =  home_dir + "hg19.fa"
        TF_RBP_ids = home_dir + "TF_RBP_human.ids"
        add_snp = True
    elif cell_type == "jurkatt":
        genome_version = "GRCh37.p13"
        species = "Human"
        ucsc_session = "http://genome.ucsc.edu/cgi-bin/hgTracks?hgS_doOtherUser=submit&hgS_otherUserName=euplotidUCSC&hgS_otherUserSessionName=primed"
        chain_file = home_dir + "hg19ToHg38.over.chain.gz"
        hsa_3D_loc = home_dir + "hES_hsa_400kb_xyz.bed"
        dna_ints = home_dir + "Jurkat_rep1and2_hg19_n3FDR.01.bedpe" #Jurkat-SMC1-petcount-filtered.bedpe,Jurkat_rep1and2_hg19_n3FDR.01.bedpe,GSE63525_CH12-LX_HiCCUPS_looplist.bed
        enh_peaks = home_dir + "Jurkat_H3K27ac_JR_peaks.bed" #Jurkat_H3K27ac_JR_peaks.bed
        prom_peaks = home_dir + "Jurkat_H3K27ac_JR_peaks.bed" #Jurkat_H3K27ac_JR_peaks.bed
        ctcf_peaks = home_dir + "20150124_CTCF_peaks.bed" #20150124_CTCF_peaks.bed
        open_peaks = home_dir + "GSE47753_CD4+_ATACseq_hg19_peaks_readcnts.bed" #GSE47753_CD4+_ATACseq_hg19_peaks_readcnts.bed
        rna_seq = home_dir + "ENCFF000MSR_jurkat_rna-seq.bed" #ENCFF000MSR_jurkat_rna-seq.bed
        encode_bed = home_dir + "hg19_encode_TFBSv2.bed"
        genome_fa =  home_dir + "hg19.fa"
        TF_RBP_ids = home_dir + "TF_RBP_human.ids"
        add_snp = True
    elif cell_type == "k562":
        genome_version = "GRCh37.p13"
        species = "Human"
        ucsc_session = "http://genome.ucsc.edu/cgi-bin/hgTracks?hgS_doOtherUser=submit&hgS_otherUserName=euplotidUCSC&hgS_otherUserSessionName=primed"
        chain_file = home_dir + "hg19ToHg38.over.chain.gz"
        hsa_3D_loc = home_dir + "hES_hsa_400kb_xyz.bed"
        dna_ints = home_dir + "GSE63525_K562_HiCCUPS_looplist.bed" #GSE63525_K562_HiCCUPS_looplist.bed,K562-RAD21-petcount-filtered.bedpe
        enh_peaks = home_dir + "ENCFF252DWA_k562_h3k27ac_hg19_narrowpeaks.bed" #ENCFF252DWA_k562_h3k27ac_hg19_narrowpeaks.bed
        prom_peaks = home_dir + "ENCFF001VCR_k562_H3K4me3_narrowpeaks.bed" #ENCFF001VCR_k562_H3K4me3_narrowpeaks.bed
        ctcf_peaks = home_dir + "ENCFF001VMZ_k562_CTCF_narrowpeak_hg19.bed" #ENCFF001VMZ_k562_CTCF_narrowpeak_hg19.bed
        open_peaks = home_dir + "ENCFF567LWS_ENCFF422QRZ_k562_DNAsePeak_hg19.bed" #ENCFF567LWS_ENCFF422QRZ_k562_DNAsePeak_hg19.bed,ENCFF001WBE_k562_DNase_narrowpeak_hg19.bed
        rna_seq = home_dir + "ENCFF000DAR_K562_rnaseq_hg19.bed" #ENCFF000DAR_K562_rnaseq_hg19.bed,ENCFF666RBZ_k562_cage_ucsc_hg19.bed
        encode_bed = home_dir + "hg19_encode_TFBSv2.bed"
        genome_fa =  home_dir + "hg19.fa"
        TF_RBP_ids = home_dir + "TF_RBP_human.ids"
        add_snp = True
    elif cell_type == "hela":
        genome_version = "GRCh37.p13"
        species = "Human"
        ucsc_session = "http://genome.ucsc.edu/cgi-bin/hgTracks?hgS_doOtherUser=submit&hgS_otherUserName=euplotidUCSC&hgS_otherUserSessionName=primed"
        chain_file = home_dir + "hg19ToHg38.over.chain.gz"
        hsa_3D_loc = home_dir + "hES_hsa_400kb_xyz.bed"
        dna_ints = home_dir + "GSE63525_HeLa_HiCCUPS_looplist.bed" #GSE63525_HeLa_HiCCUPS_looplist.bed
        enh_peaks = home_dir + "ENCFF252DWA_k562_h3k27ac_hg19_narrowpeaks.bed" #ENCFF252DWA_k562_h3k27ac_hg19_narrowpeaks.bed
        prom_peaks = home_dir + "ENCFF001VCR_k562_H3K4me3_narrowpeaks.bed" #ENCFF001VCR_k562_H3K4me3_narrowpeaks.bed
        ctcf_peaks = home_dir + "ENCFF001VMZ_k562_CTCF_narrowpeak_hg19.bed" #ENCFF001VMZ_k562_CTCF_narrowpeak_hg19.bed
        open_peaks = home_dir + "ENCFF567LWS_ENCFF422QRZ_k562_DNAsePeak_hg19.bed" #ENCFF567LWS_ENCFF422QRZ_k562_DNAsePeak_hg19.bed
        rna_seq = home_dir + "ENCFF000DAR_K562_rnaseq_hg19.bed" #ENCFF000DAR_K562_rnaseq_hg19.bed,ENCFF666RBZ_k562_cage_ucsc_hg19.bed
        encode_bed = home_dir + "hg19_encode_TFBSv2.bed"
        genome_fa =  home_dir + "hg19.fa"
        TF_RBP_ids = home_dir + "TF_RBP_human.ids"
        add_snp = True
    elif cell_type == "mesc":
        genome_version = "mm9"
        species = "Mouse"
        ucsc_session = "http://genome.ucsc.edu/cgi-bin/hgTracks?hgS_doOtherUser=submit&hgS_otherUserName=euplotidUCSC&hgS_otherUserSessionName=mesc"
        chain_file = home_dir + "mm9ToMm10.over.chain.gz"
        hsa_3D_loc = home_dir + "mES_hsa_400kb_xyz.bed"
        dna_ints = home_dir + "mES-SMC1-merged-petcount-filtered.bedpe" #mES-SMC1-merged-petcount-filtered.bedpe,Dowen_mesc_SMC1_chia_ints.bedpe,mESC_SMC_HiChip.bedpe
        enh_peaks = home_dir + "GSM2108724_chipseq_2cell_k27ac.bedGraph" #GSM2108724_chipseq_2cell_k27ac.bedGraph
        prom_peaks = home_dir + "GSM2108725_chipseq_2cell_k27me3.bedGraph" #GSM2108725_chipseq_2cell_k27me3.bedGraph
        ctcf_peaks = home_dir + "GSM747534_filt5_merged.bed" #ENCFF854IVF_mm9_mESC_CTCF.bed,GSM747534_filt5_merged.bed
        open_peaks = home_dir + "GSE66581_2cell_peaks_atac.bed" #GSE66581_2cell_peaks_atac.bed
        rna_seq = home_dir + "GSE66582_2cell_rnaseq_mm9.bed" #GSE66582_2cell_rnaseq_mm9.bed
        encode_bed = home_dir + "full_enc_mm9_tfbs.bed"
        genome_fa =  home_dir + "mm9.fa"
        TF_RBP_ids = home_dir + "TF_RBPs_mouse.ids"
        add_snp = False
    else:
        raise NameError("Cant find that cell type, try: primed, naive, neuron, NPC, GM12878, Jurkatt, K562, HeLa, mESC")

    #get list of included chromosomes from genome fasta
    homo_gen = pysam.FastaFile(genome_fa)
    #keep only canonical names chr1-22,X,Y,M
    chr_names = homo_gen.references
    chr_names = [chr for chr in chr_names if "_" not in chr]
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
    dna_int_graph = label_gene_nodes(dna_int_graph, rna_seq)
    #add TFBS from Chip-Seq under nodes
    dna_int_graph = add_TFBS_graph(dna_int_graph, encode_bed)
    #get target nodes and target node names
    target_nodes = get_target_nodes(dna_int_graph, args, cell_type)
    target_communities = set()
    for node in target_nodes:
        target_communities |=  set((dna_int_graph.node[node]["name"]).split(","))
    #Remove target communities that already have JSONs in output file
    print("Number of Target Genes: " + str(len(target_communities)))
    print("Target Genes to build Insulated Neighborhoods: " + str(",".join(target_communities)))
    #Sort nodes by FPKM for louvian
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
    num_res_steps = 20

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
            gene_nodes = [node for node in nodes if "gene" in dna_int_graph.node[node]]
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
                in_name_all_genes = ",".join([dna_int_graph.node[target_gene]["name"] for target_gene in target_overlap])
                for target_gene in target_overlap:
                    dna_int_graph.node[target_gene]["in_name"] = in_name_all_genes
                    dna_int_graph.node[target_gene]["in_min"] = min_comm
                    dna_int_graph.node[target_gene]["in_max"] = max_comm
                cur_comm_limits.write(ref_node_arr[0] + "\t" + str(min_comm) + "\t" + str(max_comm+1) + "\t" + target_overlap[0] + "\n")
        cur_comm_limits.close()
        os.system("sort-bed " + comm_out_name + " > t.bed; mv t.bed " + comm_out_name)
        #count how many of target communities fell withing CTCF - CTCF loop
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
    dna_int_graph.graph["top_genes"] = ""
    json_draft_name = (out_dir + run_name)[0:80]+ "_global_graph.json"
    json_comm = open(json_draft_name,"w+")
    json_out_graph = json_graph.node_link_data(draft_genome_out)
    json_out = json.dumps(json_out_graph)
    json_comm.write(json_out)
    json_comm.close()

    #Remove already built INs and repetitive nodes
    final_in_names = set()
    final_in_nodes = set()
    for target_node in target2partition.keys():
        in_name = dna_int_graph.node[target_node]["in_name"]
        if in_name in final_in_names:
            continue
        else:
            final_in_names.add(in_name)
            final_in_nodes.add(target_node)
    print("Annotating " + str(len(final_in_names)) + " Neighborhoods")
    print("Annotating the following Insulated Neighborhoods: " + "|".join(final_in_names))
    #Assign each open chromatin accessibility peak to the closest cohesin peak
    open_peaks_file="open_peaks_nodes.bed"
    os.system("closest-features --closest --dist " + open_peaks + " all_nodes_sorted.bed > " + open_peaks_file)
    
    num_cores = multiprocessing.cpu_count()

    results = Parallel(n_jobs=num_cores,temp_folder=out_dir) (delayed(build_IN_levels) (target_node_name, dendogram_com, dna_int_graph, genome_fa, chain_file, chain_file2, out_dir, run_name, json_draft_name,
                                                                    TF_RBP_ids, open_peaks_file, picked_tissue, ctcf2bounds, hsa_3D_loc, ucsc_session) for target_node_name in final_in_nodes)

    #Add check if succeded?
    #Annotate draft genome    
    json_comm = open(json_draft_name,"r")
    json_in = json.load(json_comm)
    json_comm.close()
    draft_genome_out_annot = json_graph.node_link_graph(json_in)

    #add colored nodes by chr
    with open("all_INs.bed","w") as all_INs:
        for node in draft_genome_out_annot.nodes():
            arr = re.split(r"[-:]",node)
            draft_genome_out_annot.node[node]["name"] = str(node) + " " + dna_int_graph.node[node]["in_name"]
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
    #dicts for number of times TF seen vs seen dimerizing, in-silico and Chip-Seq
    tf_overlap = dict()
    tf_diff = dict()
    tf_seen = set()
    #Encode peaks
    tf_overlap_encode = dict()
    tf_diff_encode = dict()
    tf_seen_encode = set()

    with open("all_INs_sorted.bed","r") as all_INs_sorted:
        for IN_node in all_INs_sorted:
            #set size of nodes according to sum of open peaks
            node_arr = IN_node.strip().split()
            node_name = (str(node_arr[0]) + ":" + str(node_arr[1]) + "-" + str(node_arr[2]))
            with open(draft_genome_out_annot.node[node_name]["json_name"]) as json_IN_file:
                json_in = json.load(json_IN_file)
                cur_IN_graph = json_graph.node_link_graph(json_in)

            #add to predicted co-occurences from sequence predictions
            for edge in cur_IN_graph.edges():
                if ("filt_tf_peak" in  cur_IN_graph.node[edge[0]]) and ("filt_tf_peak" in  cur_IN_graph.node[edge[1]]): 
                    l_node_names = set(cur_IN_graph.node[edge[0]]["filt_tf_peak"].strip().split("|"))
                    r_node_names = set(cur_IN_graph.node[edge[1]]["filt_tf_peak"].strip().split("|"))
                    over_tf_names = l_node_names.intersection(r_node_names)
                    union_tf_names = l_node_names.union(r_node_names)
                    tf_seen |= union_tf_names
                    diff_tf_names = union_tf_names.difference(over_tf_names)
                    for d_tf in diff_tf_names:
                        if d_tf not in tf_diff:
                            tf_diff[d_tf] = 0.0
                        tf_diff[d_tf] += 1.0
                    for o_tf in over_tf_names:
                        if o_tf not in tf_overlap:
                            tf_overlap[o_tf] = 0.0
                        tf_overlap[o_tf] += 1.0
                #add encode peak overlaps
                if ("tfbs_enc_all" in  cur_IN_graph.node[edge[0]]) and ("tfbs_enc_all" in  cur_IN_graph.node[edge[1]]): 
                    l_node_names = set(cur_IN_graph.node[edge[0]]["tfbs_enc_all"].strip().split("|"))
                    r_node_names = set(cur_IN_graph.node[edge[1]]["tfbs_enc_all"].strip().split("|"))
                    over_tf_names = l_node_names.intersection(r_node_names)
                    union_tf_names = l_node_names.union(r_node_names)
                    tf_seen_encode |= union_tf_names
                    diff_tf_names = union_tf_names.difference(over_tf_names)
                    for d_tf in diff_tf_names:
                        if d_tf not in tf_diff_encode:
                            tf_diff_encode[d_tf] = 0.0
                        tf_diff_encode[d_tf] += 1.0
                    for o_tf in over_tf_names:
                        if o_tf not in tf_overlap_encode:
                            tf_overlap_encode[o_tf] = 0.0
                        tf_overlap_encode[o_tf] += 1.0
            #sum size of open regions
            sum_open = 0
            for node in cur_IN_graph.nodes(data=True):
                if "open_size" in node[1]:
                    sum_open += node[1]["open_size"]
            draft_genome_out_annot.node[node_name]["open_size"] = sum_open
            draft_genome_out_annot.node[node_name]["size"] = 20
            draft_genome_out_annot.node[node_name]["color"] = color_chrom[node_arr[0]]
            sorted_ins_list.append(node_name)

    tf_dimer_perc = dict()
    for tf in tf_seen:
        n_over = 0
        n_diff = 0
        if tf in tf_overlap:
            n_over = tf_overlap[tf]
        if tf in tf_diff:
            n_diff = tf_diff[tf]
        tf_dimer_perc[tf] = int(n_over/(n_over+n_diff)*100)

    #Sort tf names by percentage seen overlapping
    tf_names = tf_dimer_perc.keys()
    tf_dimer_int = tf_dimer_perc.values()
    ix_dimer = np.array(tf_dimer_int).argsort()[::-1]
    tf_dimer_names = [tf_names[ix] for ix in ix_dimer]
    tf_dimer_perc = [str(round(tf_dimer_int[ix],4)) for ix in ix_dimer]
    draft_genome_out_annot.graph["tf_dimer_names_top100"] = ",".join(tf_dimer_names[0:min(100,len(tf_dimer_names))])
    draft_genome_out_annot.graph["tf_dimer_perc_top100"] = ",".join(tf_dimer_perc[0:min(100,len(tf_dimer_perc))])

    #Repeat for encode TFBS
    tf_dimer_perc_encode = dict()
    for tf in tf_seen_encode:
        n_over = 0
        n_diff = 0
        if tf in tf_overlap_encode:
            n_over = tf_overlap_encode[tf]
        if tf in tf_diff_encode:
            n_diff = tf_diff_encode[tf]
        tf_dimer_perc_encode[tf] = int(n_over/(n_over+n_diff)*100)

    #Sort tf names by percentage seen overlapping
    tf_names_encode = tf_dimer_perc_encode.keys()
    tf_dimer_int_encode = tf_dimer_perc_encode.values()
    ix_dimer_encode = np.array(tf_dimer_int_encode).argsort()[::-1]
    tf_dimer_names_encode = [tf_names_encode[ix] for ix in ix_dimer_encode]
    tf_dimer_perc_encode = [str(round(tf_dimer_int_encode[ix],4)) for ix in ix_dimer_encode]
    draft_genome_out_annot.graph["tf_dimer_names_top100_encode"] = ",".join(tf_dimer_names_encode[0:min(100,len(tf_dimer_names_encode))])
    draft_genome_out_annot.graph["tf_dimer_perc_top100_encode"] = ",".join(tf_dimer_perc_encode[0:min(100,len(tf_dimer_perc_encode))])

    #Add chromosome edges connecting INs 
    for i in range(0,len(sorted_ins_list)-2):
        IN1_arr = re.split(r"[-:]",sorted_ins_list[i])
        IN2_arr = re.split(r"[-:]",sorted_ins_list[i+1])
        if IN1_arr[0] == IN2_arr[0]:
            draft_genome_out_annot.add_edge(sorted_ins_list[i],sorted_ins_list[i+1])

    #Save in JSON output
    json_out_graph = json_graph.node_link_data(draft_genome_out_annot)
    json_out = json.dumps(json_out_graph)
    json_draft_annot_name = (out_dir + run_name)[0:80]+ "_global_graph_annotated.json"
    json_draft_annot = open(json_draft_annot_name,"w+")
    with open(json_draft_annot_name,"w+") as json_draft_annot:
        json_draft_annot.write(json_out)

    os.remove("all_nodes_sorted.bed")
    os.remove("open_peaks_nodes.bed")
    os.remove("ctcf_loops.bed")
    os.remove("all_INs_sorted.bed")
