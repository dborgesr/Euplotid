
# Build and annotate INs
Shows how to use the pipeline that first builds Insulated Neighborhoods using a Louvian Graph based partitioning algorithm and then annotates all Cis-regulatory elements using a Neural Network. 


```python
!python /root/Euplotid/pipelines/buildAnnotatteINs.py --help
```

    usage: buildAnnotatteINs.py [-h] -i INPUT_DIRECTORY -o OUTPUT_DIRECTORY
                                [-t INSULATED_NEIGH] -c DNA_INT_CELL_TYPE -q
                                EQTL_TISSUE [-g TOP_GENE_FPKM] -n NAME_OUTPUT
                                [-m MARKER_REGION]
    
    Pipeline to output graphical genomic model from 2D data
    
    optional arguments:
      -h, --help            show this help message and exit
      -i INPUT_DIRECTORY, --input_directory INPUT_DIRECTORY
                            Where all the input data lives
      -o OUTPUT_DIRECTORY, --output_directory OUTPUT_DIRECTORY
                            Output directory of graphs
      -t INSULATED_NEIGH, --insulated_neigh INSULATED_NEIGH
                            Comma separated list of Insulated Neighborhoods to
                            build
      -c DNA_INT_CELL_TYPE, --dna_int_cell_type DNA_INT_CELL_TYPE
                            Cell type to use in building Insulated Neighborhoods
      -q EQTL_TISSUE, --eQTL_tissue EQTL_TISSUE
                            What tissue from GTeX to do eQTL analysis w/.
                            Available tissues:
                            Adipose_Subcutaneous,Adipose_Visceral_Omentum,
                            Adrenal_Gland, Artery_Aorta, Artery_Coronary,
                            Artery_Tibial, Brain_Anterior_cingulate_cortex_BA24,
                            Brain_Caudate_basal_ganglia,
                            Brain_Cerebellar_Hemisphere, Brain_Cerebellum,
                            Brain_Cortex, Brain_Frontal_Cortex_BA9,
                            Brain_Hippocampus, Brain_Hypothalamus,
                            Brain_Nucleus_accumbens_basal_ganglia,
                            Brain_Putamen_basal_ganglia,Breast_Mammary_Tissue,
                            Cells_EBV-transformed_lymphocytes,
                            Cells_Transformed_fibroblasts, Colon_Sigmoid,
                            Colon_Transverse, Esophagus_Gastroesophageal_Junction,
                            Esophagus_Mucosa, Esophagus_Muscularis,
                            Heart_Atrial_Appendage, Heart_Left_Ventricle, Liver,
                            Lung, Muscle_Skeletal, Nerve_Tibial, Ovary, Pancreas,
                            Pituitary, Prostate, Skin_Not_Sun_Exposed_Suprapubic,
                            Skin_Sun_Exposed_Lower_leg,
                            Small_Intestine_Terminal_Ileum, Spleen, Stomach,
                            Testis, Thyroid, Uterus, Vagina, Whole_Blood
      -g TOP_GENE_FPKM, --top_gene_fpkm TOP_GENE_FPKM
                            Instead of defining Insulated Neighborhoods take genes
                            expressed above FPKM cutoff
      -n NAME_OUTPUT, --name_output NAME_OUTPUT
                            Output prefix for all output files
      -m MARKER_REGION, --marker_region MARKER_REGION
                            Build all Insulated Neighborhoods of genes which fall
                            within this chromosomal region (ex:
                            chr1:54996039-55996039)



```python
!luarocks install trepl
```

    
    Error: No results matching query were found.



```python
!python /root/Euplotid/src/buildAnnotateINs.py \
    -i /input_dir/ -o /output_dir/ \
    -c primed -q Adipose_Visceral_Omentum \
    -n example -t IRX3,FTO,CHD9
```

    Number of Target Genes: 6
    Target Genes to build Insulated Neighborhoods: RPGRIP1L,JB149426,FTO,Y_RNA,IRX3,CHD9
    Checking resolution parameter against CTCF-CTCF loop
    Annotating 4 Neighborhoods
    Annotating the following Insulated Neighborhoods: CHD9|RPGRIP1L,FTO|CHD9,CHD9,CHD9,Y_RNA,CHD9|IRX3,FTO,JB149426,FTO
    Annotating Insulated Neighborhood: IRX3,FTO,JB149426,FTO
    Annotating Open Regions
    Annotating Insulated Neighborhood: CHD9
    Annotating Open Regions
    Annotating Insulated Neighborhood: CHD9,CHD9,CHD9,Y_RNA,CHD9
    Annotating Open Regions
    Annotating Insulated Neighborhood: RPGRIP1L,FTO
    Annotating Open Regions
    Annotating SNPs/CNVs
    querying 1-4...done.
    Annotating SNPs/CNVs
    querying 1-1...done.
    Annotating SNPs/CNVs
    querying 1-5...done.
    Annotating SNPs/CNVs
    querying 1-2...done.
    Fetching 12 variant(s) . . .
    Fetching 17 variant(s) . . .
    Fetching 52 variant(s) . . .
    Fetching 18 variant(s) . . .
    No results to return
    Fetching 29 variant(s) . . .
    No results to return
    Fetching 34 variant(s) . . .
    No results to return
    Fetching 49 variant(s) . . .
    No results to return
    Deepbind SNP prediction of: chr16:53979176-53980851
    No results to return
    Deepbind SNP prediction of: chr16:53241256-53242869
    No results to return
    Deepbind SNP prediction of: chr16:53254174-53255315
    No results to return
    Fetching 54 variant(s) . . .
    No results to return
    Fetching 77 variant(s) . . .
    No results to return
    Fetching 23 variant(s) . . .
    No results to return
    Fetching 83 variant(s) . . .
    No results to return
    Fetching 608 variant(s) . . .
    No results to return
    Fetching 58 variant(s) . . .
    No results to return
    Fetching 88 variant(s) . . .
    No results to return
    Fetching 2393 variant(s) . . .
    Basset SNP prediction of: chr16:53979176-53980851
    /root/Basset/src/basset_sad.py -f /input_dir/hg19.fa -l 600 -o /output_dir//primed/example//chr16_53979176_53980851_basset -t /root/Basset/tutorials/sad_eg/sample_beds.txt /input_dir/pretrained_model.th /output_dir//primed/example//chr16_53979176_53980851_basset/open_variants.vcf
    No results to return
    Fetching 25 variant(s) . . .
    WARNING: Skipping chr16:g.53932919_53932920del - neither allele matches reference genome: TAA vs AAA
    WARNING: Skipping chr16:g.53932930_53932932del - neither allele matches reference genome: CTCT vs TCTT
    WARNING: Skipping chr16:g.53933428del - neither allele matches reference genome: TC vs CC
    /root/Basset/src/basset_predict.lua -rc  /input_dir/pretrained_model.th /output_dir//primed/example//chr16_53979176_53980851_basset/model_in.h5 /output_dir//primed/example//chr16_53979176_53980851_basset/model_out.txt
    /opt/conda/bin/lua: /opt/conda/share/lua/5.2/luarocks/loader.lua:117: error loading module 'treplutils' from file '/opt/conda/bin/lua':
    	/opt/conda/bin/lua:1: unexpected symbol near char(127)
    stack traceback:
    	[C]: in function 'a_loader'
    	/opt/conda/share/lua/5.2/luarocks/loader.lua:117: in function </opt/conda/share/lua/5.2/luarocks/loader.lua:114>
    	(...tail calls...)
    	[C]: in function 'require'
    	/opt/conda/share/lua/5.2/trepl/init.lua:40: in main chunk
    	[C]: in function 'require'
    	/opt/conda/lib/luarocks/rocks/trepl/scm-1/bin/th:104: in main chunk
    	[C]: in ?
    Traceback (most recent call last):
      File "/root/Basset/src/basset_sad.py", line 222, in <module>
        main()
      File "/root/Basset/src/basset_sad.py", line 92, in main
        for line in open(options.model_hdf5_file):
    FileNotFoundError: [Errno 2] No such file or directory: '/output_dir//primed/example//chr16_53979176_53980851_basset/model_out.txt'
    multiprocessing.pool.RemoteTraceback: 
    """
    Traceback (most recent call last):
      File "/opt/conda/lib/python3.5/site-packages/joblib/_parallel_backends.py", line 350, in __call__
        return self.func(*args, **kwargs)
      File "/opt/conda/lib/python3.5/site-packages/joblib/parallel.py", line 131, in __call__
        return [func(*args, **kwargs) for func, args, kwargs in self.items]
      File "/opt/conda/lib/python3.5/site-packages/joblib/parallel.py", line 131, in <listcomp>
        return [func(*args, **kwargs) for func, args, kwargs in self.items]
      File "/root/Euplotid/src/buildAnnotateINs.py", line 781, in build_IN_levels
        comm_var_out, comm_nodes_out, comm_edges_out, comm_edges_out_bed, ucsc_session_in,out_dir)
      File "/root/Euplotid/src/buildAnnotateINs.py", line 673, in annotate_IN
        source_community_graph = add_variants_predict(source_community_graph, homo_gen, chain_file, TF_RBP_ids, picked_tissue, target_node_name,out_dir)
      File "/root/Euplotid/src/buildAnnotateINs.py", line 302, in add_variants_predict
        G = add_basset_sad_sat(G, homo_gen, target_node_name,out_dir)
      File "/root/Euplotid/src/buildAnnotateINs.py", line 406, in add_basset_sad_sat
        sad_table = pd.read_table(basset_out_name + "/sad_table.txt", delim_whitespace=True, header = 0)
      File "/opt/conda/lib/python3.5/site-packages/pandas/io/parsers.py", line 655, in parser_f
        return _read(filepath_or_buffer, kwds)
      File "/opt/conda/lib/python3.5/site-packages/pandas/io/parsers.py", line 405, in _read
        parser = TextFileReader(filepath_or_buffer, **kwds)
      File "/opt/conda/lib/python3.5/site-packages/pandas/io/parsers.py", line 762, in __init__
        self._make_engine(self.engine)
      File "/opt/conda/lib/python3.5/site-packages/pandas/io/parsers.py", line 966, in _make_engine
        self._engine = CParserWrapper(self.f, **self.options)
      File "/opt/conda/lib/python3.5/site-packages/pandas/io/parsers.py", line 1582, in __init__
        self._reader = parsers.TextReader(src, **kwds)
      File "pandas/_libs/parsers.pyx", line 394, in pandas._libs.parsers.TextReader.__cinit__ (pandas/_libs/parsers.c:4209)
      File "pandas/_libs/parsers.pyx", line 710, in pandas._libs.parsers.TextReader._setup_parser_source (pandas/_libs/parsers.c:8873)
    FileNotFoundError: File b'/output_dir//primed/example//chr16_53979176_53980851_basset/sad_table.txt' does not exist
    
    During handling of the above exception, another exception occurred:
    
    Traceback (most recent call last):
      File "/opt/conda/lib/python3.5/multiprocessing/pool.py", line 119, in worker
        result = (True, func(*args, **kwds))
      File "/opt/conda/lib/python3.5/site-packages/joblib/_parallel_backends.py", line 359, in __call__
        raise TransportableException(text, e_type)
    joblib.my_exceptions.TransportableException: TransportableException
    ___________________________________________________________________________
    FileNotFoundError                                  Wed May 31 03:20:16 2017
    PID: 3714                               Python 3.5.3: /opt/conda/bin/python
    ...........................................................................
    /opt/conda/lib/python3.5/site-packages/joblib/parallel.py in __call__(self=<joblib.parallel.BatchedCalls object>)
        126     def __init__(self, iterator_slice):
        127         self.items = list(iterator_slice)
        128         self._size = len(self.items)
        129 
        130     def __call__(self):
    --> 131         return [func(*args, **kwargs) for func, args, kwargs in self.items]
            self.items = [(<function build_IN_levels>, ('chr16:53979176-53980851', [{'chr10:100002276-100003364': 108, 'chr10:100021906-100023080': 108, 'chr10:100067964-100070284': 108, 'chr10:100091165-100092468': 108, 'chr10:100174343-100175470': 108, 'chr10:100183011-100184125': 1475, 'chr10:100185094-100186268': 1475, 'chr10:100205605-100207134': 1475, 'chr10:100226542-100228626': 1475, 'chr10:100423690-100425020': 1475, ...}, {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7, 8: 8, 9: 9, ...}, {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7, 8: 8, 9: 9, ...}, {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7, 8: 8, 9: 9, ...}], <networkx.classes.graph.Graph object>, '/input_dir/hg19.fa', '/input_dir/hg19ToHg38.over.chain.gz', '/input_dir/hg38ToHg19.over.chain.gz', '/output_dir//primed/example/', 'example', '/output_dir//primed/example/example_global_graph.json', '/input_dir/TF_RBP_human.ids', 'open_peaks_nodes.bed', 'Adipose_Visceral_Omentum', {'chr16:53088027-53090207': ('chr16:53088027-53090207', 'chr16:86527974-86528787'), 'chr16:53132274-53134474': ('chr16:53088027-53090207', 'chr16:86527974-86528787'), 'chr16:53193391-53194589': ('chr16:53088027-53090207', 'chr16:53193391-53194589'), 'chr16:53241256-53242869': ('chr16:53055849-53057134', 'chr16:52101824-52104715'), 'chr16:53254174-53255315': ('chr16:53088027-53090207', 'chr16:53254174-53255315'), 'chr16:53776635-53777704': ('chr5:174007939-174009432', 'chr5:130752381-130753830'), 'chr16:53925426-53926800': ('chr16:54315183-54322626', 'chr16:53925426-53926800'), 'chr16:53979176-53980851': ('chr16:54315183-54322626', 'chr16:53925426-53926800'), 'chr16:54315183-54322626': ('chr16:53088027-53090207', 'chr16:86527974-86528787')}, '/input_dir/hES_hsa_400kb_xyz.bed', 'http://genome.ucsc.edu/cgi-bin/hgTracks?hgS_doOt...Name=euplotidUCSC&hgS_otherUserSessionName=primed'), {})]
        132 
        133     def __len__(self):
        134         return self._size
        135 
    
    ...........................................................................
    /opt/conda/lib/python3.5/site-packages/joblib/parallel.py in <listcomp>(.0=<list_iterator object>)
        126     def __init__(self, iterator_slice):
        127         self.items = list(iterator_slice)
        128         self._size = len(self.items)
        129 
        130     def __call__(self):
    --> 131         return [func(*args, **kwargs) for func, args, kwargs in self.items]
            func = <function build_IN_levels>
            args = ('chr16:53979176-53980851', [{'chr10:100002276-100003364': 108, 'chr10:100021906-100023080': 108, 'chr10:100067964-100070284': 108, 'chr10:100091165-100092468': 108, 'chr10:100174343-100175470': 108, 'chr10:100183011-100184125': 1475, 'chr10:100185094-100186268': 1475, 'chr10:100205605-100207134': 1475, 'chr10:100226542-100228626': 1475, 'chr10:100423690-100425020': 1475, ...}, {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7, 8: 8, 9: 9, ...}, {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7, 8: 8, 9: 9, ...}, {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7, 8: 8, 9: 9, ...}], <networkx.classes.graph.Graph object>, '/input_dir/hg19.fa', '/input_dir/hg19ToHg38.over.chain.gz', '/input_dir/hg38ToHg19.over.chain.gz', '/output_dir//primed/example/', 'example', '/output_dir//primed/example/example_global_graph.json', '/input_dir/TF_RBP_human.ids', 'open_peaks_nodes.bed', 'Adipose_Visceral_Omentum', {'chr16:53088027-53090207': ('chr16:53088027-53090207', 'chr16:86527974-86528787'), 'chr16:53132274-53134474': ('chr16:53088027-53090207', 'chr16:86527974-86528787'), 'chr16:53193391-53194589': ('chr16:53088027-53090207', 'chr16:53193391-53194589'), 'chr16:53241256-53242869': ('chr16:53055849-53057134', 'chr16:52101824-52104715'), 'chr16:53254174-53255315': ('chr16:53088027-53090207', 'chr16:53254174-53255315'), 'chr16:53776635-53777704': ('chr5:174007939-174009432', 'chr5:130752381-130753830'), 'chr16:53925426-53926800': ('chr16:54315183-54322626', 'chr16:53925426-53926800'), 'chr16:53979176-53980851': ('chr16:54315183-54322626', 'chr16:53925426-53926800'), 'chr16:54315183-54322626': ('chr16:53088027-53090207', 'chr16:86527974-86528787')}, '/input_dir/hES_hsa_400kb_xyz.bed', 'http://genome.ucsc.edu/cgi-bin/hgTracks?hgS_doOt...Name=euplotidUCSC&hgS_otherUserSessionName=primed')
            kwargs = {}
        132 
        133     def __len__(self):
        134         return self._size
        135 
    
    ...........................................................................
    /root/Euplotid/src/buildAnnotateINs.py in build_IN_levels(target_node_name='chr16:53979176-53980851', dendogram_com=[{'chr10:100002276-100003364': 108, 'chr10:100021906-100023080': 108, 'chr10:100067964-100070284': 108, 'chr10:100091165-100092468': 108, 'chr10:100174343-100175470': 108, 'chr10:100183011-100184125': 1528, 'chr10:100185094-100186268': 1528, 'chr10:100205605-100207134': 7231, 'chr10:100226542-100228626': 7231, 'chr10:100423690-100425020': 1528, ...}, {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7, 8: 8, 9: 9, ...}, {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7, 8: 8, 9: 9, ...}], dna_int_graph=<networkx.classes.graph.Graph object>, genome_fa='/input_dir/hg19.fa', chain_file='/input_dir/hg19ToHg38.over.chain.gz', chain_file2='/input_dir/hg38ToHg19.over.chain.gz', out_dir='/output_dir//primed/example/', run_name='example', json_draft_name='/output_dir//primed/example/example_global_graph.json', TF_RBP_ids='/input_dir/TF_RBP_human.ids', open_peaks_file='open_peaks_nodes.bed', picked_tissue='Adipose_Visceral_Omentum', ctcf2bounds={'chr16:53088027-53090207': ('chr16:53088027-53090207', 'chr16:86527974-86528787'), 'chr16:53132274-53134474': ('chr16:53088027-53090207', 'chr16:86527974-86528787'), 'chr16:53193391-53194589': ('chr16:53088027-53090207', 'chr16:53193391-53194589'), 'chr16:53241256-53242869': ('chr16:53055849-53057134', 'chr16:52101824-52104715'), 'chr16:53254174-53255315': ('chr16:53088027-53090207', 'chr16:53254174-53255315'), 'chr16:53776635-53777704': ('chr5:174007939-174009432', 'chr5:130752381-130753830'), 'chr16:53925426-53926800': ('chr16:54315183-54322626', 'chr16:53925426-53926800'), 'chr16:53979176-53980851': ('chr16:54315183-54322626', 'chr16:53925426-53926800'), 'chr16:54315183-54322626': ('chr16:53088027-53090207', 'chr16:86527974-86528787')}, hsa_3D_loc='/input_dir/hES_hsa_400kb_xyz.bed', ucsc_session='http://genome.ucsc.edu/cgi-bin/hgTracks?hgS_doOt...Name=euplotidUCSC&hgS_otherUserSessionName=primed')
        776         add_snp = True
        777     else:
        778         add_snp = False
        779     source_community_graph = annotate_IN(target_node_name, dendogram_com, level, dna_int_graph, add_open, add_snp,
        780                                                   homo_gen, chain_file, chain_file2, TF_RBP_ids, open_peaks_file, picked_tissue, ctcf2bounds, comm_open_out,
    --> 781                                                   comm_var_out, comm_nodes_out, comm_edges_out, comm_edges_out_bed, ucsc_session_in,out_dir)
            comm_var_out = <_io.TextIOWrapper name='/output_dir//primed/exa...y_variants_lvl_1.bed' mode='a+' encoding='UTF-8'>
            comm_nodes_out = <_io.TextIOWrapper name='/output_dir//primed/exa...nity_nodes_lvl_1.bed' mode='a+' encoding='UTF-8'>
            comm_edges_out = <_io.TextIOWrapper name='/output_dir//primed/exa...ty_edges_lvl_1.washu' mode='a+' encoding='UTF-8'>
            comm_edges_out_bed = <_io.TextIOWrapper name='/output_dir//primed/exa...nity_edges_lvl_1.bed' mode='a+' encoding='UTF-8'>
            ucsc_session_in = 'http://genome.ucsc.edu/cgi-bin/hgTracks?hgS_doOt...ssionName=primed&position=chr16:54668601-53806737'
            out_dir = '/output_dir//primed/example/'
        782     #level 1 
        783     if len(dendogram_com) > 1:
        784         level = 1
        785         add_open = False
    
    ...........................................................................
    /root/Euplotid/src/buildAnnotateINs.py in annotate_IN(target_node_name='chr16:53979176-53980851', dendogram_com=[{'chr10:100002276-100003364': 108, 'chr10:100021906-100023080': 108, 'chr10:100067964-100070284': 108, 'chr10:100091165-100092468': 108, 'chr10:100174343-100175470': 108, 'chr10:100183011-100184125': 1528, 'chr10:100185094-100186268': 1528, 'chr10:100205605-100207134': 7231, 'chr10:100226542-100228626': 7231, 'chr10:100423690-100425020': 1528, ...}, {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7, 8: 8, 9: 9, ...}, {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7, 8: 8, 9: 9, ...}], level=0, dna_int_graph=<networkx.classes.graph.Graph object>, add_open=True, add_snp=True, homo_gen=<pysam.libcfaidx.FastaFile object>, chain_file='/input_dir/hg19ToHg38.over.chain.gz', chain_file2='/input_dir/hg38ToHg19.over.chain.gz', TF_RBP_ids='/input_dir/TF_RBP_human.ids', open_peaks_file='open_peaks_nodes.bed', picked_tissue='Adipose_Visceral_Omentum', ctcf2bounds={'chr16:53088027-53090207': ('chr16:53088027-53090207', 'chr16:86527974-86528787'), 'chr16:53132274-53134474': ('chr16:53088027-53090207', 'chr16:86527974-86528787'), 'chr16:53193391-53194589': ('chr16:53088027-53090207', 'chr16:53193391-53194589'), 'chr16:53241256-53242869': ('chr16:53055849-53057134', 'chr16:52101824-52104715'), 'chr16:53254174-53255315': ('chr16:53088027-53090207', 'chr16:53254174-53255315'), 'chr16:53776635-53777704': ('chr5:174007939-174009432', 'chr5:130752381-130753830'), 'chr16:53925426-53926800': ('chr16:54315183-54322626', 'chr16:53925426-53926800'), 'chr16:53979176-53980851': ('chr16:54315183-54322626', 'chr16:53925426-53926800'), 'chr16:54315183-54322626': ('chr16:53088027-53090207', 'chr16:86527974-86528787')}, comm_open_out=<_io.TextIOWrapper name='/output_dir//primed/exa...en_regions_lvl_1.bed' mode='a+' encoding='UTF-8'>, comm_var_out=<_io.TextIOWrapper name='/output_dir//primed/exa...y_variants_lvl_1.bed' mode='a+' encoding='UTF-8'>, comm_nodes_out=<_io.TextIOWrapper name='/output_dir//primed/exa...nity_nodes_lvl_1.bed' mode='a+' encoding='UTF-8'>, comm_edges_out=<_io.TextIOWrapper name='/output_dir//primed/exa...ty_edges_lvl_1.washu' mode='a+' encoding='UTF-8'>, comm_edges_out_bed=<_io.TextIOWrapper name='/output_dir//primed/exa...nity_edges_lvl_1.bed' mode='a+' encoding='UTF-8'>, ucsc_session_in='http://genome.ucsc.edu/cgi-bin/hgTracks?hgS_doOt...ssionName=primed&position=chr16:54668601-53806737', out_dir='/output_dir//primed/example/')
        668         print("Annotating Open Regions")
        669         source_community_graph = add_openRegions_predict(source_community_graph, homo_gen, chain_file, chain_file2, TF_RBP_ids, open_peaks_file, target_node_name, level, dict())
        670         sys.stdout.flush()
        671         if add_snp:
        672             print("Annotating SNPs/CNVs")
    --> 673             source_community_graph = add_variants_predict(source_community_graph, homo_gen, chain_file, TF_RBP_ids, picked_tissue, target_node_name,out_dir)
            source_community_graph = <networkx.classes.graph.Graph object>
            homo_gen = <pysam.libcfaidx.FastaFile object>
            chain_file = '/input_dir/hg19ToHg38.over.chain.gz'
            TF_RBP_ids = '/input_dir/TF_RBP_human.ids'
            picked_tissue = 'Adipose_Visceral_Omentum'
            target_node_name = 'chr16:53979176-53980851'
            out_dir = '/output_dir//primed/example/'
        674             sys.stdout.flush()
        675     
        676     #Color Insulated Neighborhood according to node attributes, make mid distance to target node
        677 
    
    ...........................................................................
    /root/Euplotid/src/buildAnnotateINs.py in add_variants_predict(G=<networkx.classes.graph.Graph object>, homo_gen=<pysam.libcfaidx.FastaFile object>, chain_file='/input_dir/hg19ToHg38.over.chain.gz', TF_RBP_ids='/input_dir/TF_RBP_human.ids', tissue_type='Adipose_Visceral_Omentum', target_node_name='chr16:53979176-53980851', out_dir='/output_dir//primed/example/')
        297     #predict effect of all variation added to graph using DeepBind and Basset
        298     if add_any_SNPs:
        299         print("Deepbind SNP prediction of: " + target_node_name)
        300         G = deepbind_predict_SNPs(G, homo_gen, TF_RBP_ids,target_node_name)
        301         print("Basset SNP prediction of: " + target_node_name)
    --> 302         G = add_basset_sad_sat(G, homo_gen, target_node_name,out_dir)
            G = <networkx.classes.graph.Graph object>
            homo_gen = <pysam.libcfaidx.FastaFile object>
            target_node_name = 'chr16:53979176-53980851'
            out_dir = '/output_dir//primed/example/'
        303     try:
        304         os.remove(basset_name_out + "open_variants.vcf")
        305     except OSError:
        306         pass
    
    ...........................................................................
    /root/Euplotid/src/buildAnnotateINs.py in add_basset_sad_sat(G=<networkx.classes.graph.Graph object>, homo_gen=<pysam.libcfaidx.FastaFile object>, target_node_name='chr16:53979176-53980851', out_dir='/output_dir//primed/example/')
        401             count += 1
        402     cmd = ("/root/Basset/src/basset_sad.py -f " + (homo_gen.filename).decode("utf-8") + " -l 600 -o " + basset_out_name + " -t " + targets_file + " " +  model_file + " " + basset_out_name + "/open_variants.vcf")
        403     print(cmd)
        404     os.system(cmd)
        405     #Pick SNP maximizing for sum of Delta SAD across all cell types in trained model
    --> 406     sad_table = pd.read_table(basset_out_name + "/sad_table.txt", delim_whitespace=True, header = 0)
            sad_table = undefined
            basset_out_name = '/output_dir//primed/example//chr16_53979176_53980851_basset'
        407     sad_table_dense = sad_table.pivot_table(index='rsid', columns='target', values='pred')
        408     serial_sad_table = pickle.dumps(sad_table_dense, protocol=0)
        409     G.graph["sad_table"] = serial_sad_table
        410     sad_table["sad_abs"] = abs(sad_table["pred"])
    
    ...........................................................................
    /opt/conda/lib/python3.5/site-packages/pandas/io/parsers.py in parser_f(filepath_or_buffer='/output_dir//primed/example//chr16_53979176_53980851_basset/sad_table.txt', sep='\t', delimiter='\t', header=0, names=None, index_col=None, usecols=None, squeeze=False, prefix=None, mangle_dupe_cols=True, dtype=None, engine='c', converters=None, true_values=None, false_values=None, skipinitialspace=False, skiprows=None, nrows=None, na_values=None, keep_default_na=True, na_filter=True, verbose=False, skip_blank_lines=True, parse_dates=False, infer_datetime_format=False, keep_date_col=False, date_parser=None, dayfirst=False, iterator=False, chunksize=None, compression='infer', thousands=None, decimal=b'.', lineterminator=None, quotechar='"', quoting=0, escapechar=None, comment=None, encoding=None, dialect=None, tupleize_cols=False, error_bad_lines=True, warn_bad_lines=True, skipfooter=0, skip_footer=0, doublequote=True, delim_whitespace=True, as_recarray=False, compact_ints=False, use_unsigned=False, low_memory=True, buffer_lines=None, memory_map=False, float_precision=None)
        650                     mangle_dupe_cols=mangle_dupe_cols,
        651                     tupleize_cols=tupleize_cols,
        652                     infer_datetime_format=infer_datetime_format,
        653                     skip_blank_lines=skip_blank_lines)
        654 
    --> 655         return _read(filepath_or_buffer, kwds)
            filepath_or_buffer = '/output_dir//primed/example//chr16_53979176_53980851_basset/sad_table.txt'
            kwds = {'as_recarray': False, 'buffer_lines': None, 'chunksize': None, 'comment': None, 'compact_ints': False, 'compression': None, 'converters': None, 'date_parser': None, 'dayfirst': False, 'decimal': b'.', ...}
        656 
        657     parser_f.__name__ = name
        658 
        659     return parser_f
    
    ...........................................................................
    /opt/conda/lib/python3.5/site-packages/pandas/io/parsers.py in _read(filepath_or_buffer='/output_dir//primed/example//chr16_53979176_53980851_basset/sad_table.txt', kwds={'as_recarray': False, 'buffer_lines': None, 'chunksize': None, 'comment': None, 'compact_ints': False, 'compression': None, 'converters': None, 'date_parser': None, 'dayfirst': False, 'decimal': b'.', ...})
        400     iterator = kwds.get('iterator', False)
        401     chunksize = _validate_integer('chunksize', kwds.get('chunksize', None), 1)
        402     nrows = _validate_integer('nrows', kwds.get('nrows', None))
        403 
        404     # Create the parser.
    --> 405     parser = TextFileReader(filepath_or_buffer, **kwds)
            parser = undefined
            filepath_or_buffer = '/output_dir//primed/example//chr16_53979176_53980851_basset/sad_table.txt'
            kwds = {'as_recarray': False, 'buffer_lines': None, 'chunksize': None, 'comment': None, 'compact_ints': False, 'compression': None, 'converters': None, 'date_parser': None, 'dayfirst': False, 'decimal': b'.', ...}
        406 
        407     if chunksize or iterator:
        408         return parser
        409 
    
    ...........................................................................
    /opt/conda/lib/python3.5/site-packages/pandas/io/parsers.py in __init__(self=<pandas.io.parsers.TextFileReader object>, f='/output_dir//primed/example//chr16_53979176_53980851_basset/sad_table.txt', engine='c', **kwds={'as_recarray': False, 'buffer_lines': None, 'chunksize': None, 'comment': None, 'compact_ints': False, 'compression': None, 'converters': None, 'date_parser': None, 'dayfirst': False, 'decimal': b'.', ...})
        757         # might mutate self.engine
        758         self.options, self.engine = self._clean_options(options, engine)
        759         if 'has_index_names' in kwds:
        760             self.options['has_index_names'] = kwds['has_index_names']
        761 
    --> 762         self._make_engine(self.engine)
            self._make_engine = <bound method TextFileReader._make_engine of <pandas.io.parsers.TextFileReader object>>
            self.engine = 'c'
        763 
        764     def close(self):
        765         self._engine.close()
        766 
    
    ...........................................................................
    /opt/conda/lib/python3.5/site-packages/pandas/io/parsers.py in _make_engine(self=<pandas.io.parsers.TextFileReader object>, engine='c')
        961             self.close()
        962             raise
        963 
        964     def _make_engine(self, engine='c'):
        965         if engine == 'c':
    --> 966             self._engine = CParserWrapper(self.f, **self.options)
            self._engine = None
            self.f = '/output_dir//primed/example//chr16_53979176_53980851_basset/sad_table.txt'
            self.options = {'as_recarray': False, 'buffer_lines': None, 'comment': None, 'compact_ints': False, 'compression': None, 'converters': {}, 'date_parser': None, 'dayfirst': False, 'decimal': b'.', 'delim_whitespace': True, ...}
        967         else:
        968             if engine == 'python':
        969                 klass = PythonParser
        970             elif engine == 'python-fwf':
    
    ...........................................................................
    /opt/conda/lib/python3.5/site-packages/pandas/io/parsers.py in __init__(self=<pandas.io.parsers.CParserWrapper object>, src='/output_dir//primed/example//chr16_53979176_53980851_basset/sad_table.txt', **kwds={'allow_leading_cols': True, 'as_recarray': False, 'buffer_lines': None, 'comment': None, 'compact_ints': False, 'compression': None, 'converters': {}, 'decimal': b'.', 'delim_whitespace': True, 'delimiter': '\t', ...})
       1577             kwds['encoding'] = 'utf-8'
       1578 
       1579         # #2442
       1580         kwds['allow_leading_cols'] = self.index_col is not False
       1581 
    -> 1582         self._reader = parsers.TextReader(src, **kwds)
            self._reader = undefined
            src = '/output_dir//primed/example//chr16_53979176_53980851_basset/sad_table.txt'
            kwds = {'allow_leading_cols': True, 'as_recarray': False, 'buffer_lines': None, 'comment': None, 'compact_ints': False, 'compression': None, 'converters': {}, 'decimal': b'.', 'delim_whitespace': True, 'delimiter': '\t', ...}
       1583 
       1584         # XXX
       1585         self.usecols, self.usecols_dtype = _validate_usecols_arg(
       1586             self._reader.usecols)
    
    ...........................................................................
    /opt/conda/lib/python3.5/site-packages/pandas/_libs/parsers.cpython-35m-x86_64-linux-gnu.so in pandas._libs.parsers.TextReader.__cinit__ (pandas/_libs/parsers.c:4209)()
    
    ...........................................................................
    /opt/conda/lib/python3.5/site-packages/pandas/_libs/parsers.cpython-35m-x86_64-linux-gnu.so in pandas._libs.parsers.TextReader._setup_parser_source (pandas/_libs/parsers.c:8873)()
    
    FileNotFoundError: File b'/output_dir//primed/example//chr16_53979176_53980851_basset/sad_table.txt' does not exist
    ___________________________________________________________________________
    """
    
    The above exception was the direct cause of the following exception:
    
    Traceback (most recent call last):
      File "/opt/conda/lib/python3.5/site-packages/joblib/parallel.py", line 699, in retrieve
        self._output.extend(job.get(timeout=self.timeout))
      File "/opt/conda/lib/python3.5/multiprocessing/pool.py", line 608, in get
        raise self._value
    joblib.my_exceptions.TransportableException: TransportableException
    ___________________________________________________________________________
    FileNotFoundError                                  Wed May 31 03:20:16 2017
    PID: 3714                               Python 3.5.3: /opt/conda/bin/python
    ...........................................................................
    /opt/conda/lib/python3.5/site-packages/joblib/parallel.py in __call__(self=<joblib.parallel.BatchedCalls object>)
        126     def __init__(self, iterator_slice):
        127         self.items = list(iterator_slice)
        128         self._size = len(self.items)
        129 
        130     def __call__(self):
    --> 131         return [func(*args, **kwargs) for func, args, kwargs in self.items]
            self.items = [(<function build_IN_levels>, ('chr16:53979176-53980851', [{'chr10:100002276-100003364': 108, 'chr10:100021906-100023080': 108, 'chr10:100067964-100070284': 108, 'chr10:100091165-100092468': 108, 'chr10:100174343-100175470': 108, 'chr10:100183011-100184125': 1475, 'chr10:100185094-100186268': 1475, 'chr10:100205605-100207134': 1475, 'chr10:100226542-100228626': 1475, 'chr10:100423690-100425020': 1475, ...}, {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7, 8: 8, 9: 9, ...}, {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7, 8: 8, 9: 9, ...}, {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7, 8: 8, 9: 9, ...}], <networkx.classes.graph.Graph object>, '/input_dir/hg19.fa', '/input_dir/hg19ToHg38.over.chain.gz', '/input_dir/hg38ToHg19.over.chain.gz', '/output_dir//primed/example/', 'example', '/output_dir//primed/example/example_global_graph.json', '/input_dir/TF_RBP_human.ids', 'open_peaks_nodes.bed', 'Adipose_Visceral_Omentum', {'chr16:53088027-53090207': ('chr16:53088027-53090207', 'chr16:86527974-86528787'), 'chr16:53132274-53134474': ('chr16:53088027-53090207', 'chr16:86527974-86528787'), 'chr16:53193391-53194589': ('chr16:53088027-53090207', 'chr16:53193391-53194589'), 'chr16:53241256-53242869': ('chr16:53055849-53057134', 'chr16:52101824-52104715'), 'chr16:53254174-53255315': ('chr16:53088027-53090207', 'chr16:53254174-53255315'), 'chr16:53776635-53777704': ('chr5:174007939-174009432', 'chr5:130752381-130753830'), 'chr16:53925426-53926800': ('chr16:54315183-54322626', 'chr16:53925426-53926800'), 'chr16:53979176-53980851': ('chr16:54315183-54322626', 'chr16:53925426-53926800'), 'chr16:54315183-54322626': ('chr16:53088027-53090207', 'chr16:86527974-86528787')}, '/input_dir/hES_hsa_400kb_xyz.bed', 'http://genome.ucsc.edu/cgi-bin/hgTracks?hgS_doOt...Name=euplotidUCSC&hgS_otherUserSessionName=primed'), {})]
        132 
        133     def __len__(self):
        134         return self._size
        135 
    
    ...........................................................................
    /opt/conda/lib/python3.5/site-packages/joblib/parallel.py in <listcomp>(.0=<list_iterator object>)
        126     def __init__(self, iterator_slice):
        127         self.items = list(iterator_slice)
        128         self._size = len(self.items)
        129 
        130     def __call__(self):
    --> 131         return [func(*args, **kwargs) for func, args, kwargs in self.items]
            func = <function build_IN_levels>
            args = ('chr16:53979176-53980851', [{'chr10:100002276-100003364': 108, 'chr10:100021906-100023080': 108, 'chr10:100067964-100070284': 108, 'chr10:100091165-100092468': 108, 'chr10:100174343-100175470': 108, 'chr10:100183011-100184125': 1475, 'chr10:100185094-100186268': 1475, 'chr10:100205605-100207134': 1475, 'chr10:100226542-100228626': 1475, 'chr10:100423690-100425020': 1475, ...}, {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7, 8: 8, 9: 9, ...}, {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7, 8: 8, 9: 9, ...}, {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7, 8: 8, 9: 9, ...}], <networkx.classes.graph.Graph object>, '/input_dir/hg19.fa', '/input_dir/hg19ToHg38.over.chain.gz', '/input_dir/hg38ToHg19.over.chain.gz', '/output_dir//primed/example/', 'example', '/output_dir//primed/example/example_global_graph.json', '/input_dir/TF_RBP_human.ids', 'open_peaks_nodes.bed', 'Adipose_Visceral_Omentum', {'chr16:53088027-53090207': ('chr16:53088027-53090207', 'chr16:86527974-86528787'), 'chr16:53132274-53134474': ('chr16:53088027-53090207', 'chr16:86527974-86528787'), 'chr16:53193391-53194589': ('chr16:53088027-53090207', 'chr16:53193391-53194589'), 'chr16:53241256-53242869': ('chr16:53055849-53057134', 'chr16:52101824-52104715'), 'chr16:53254174-53255315': ('chr16:53088027-53090207', 'chr16:53254174-53255315'), 'chr16:53776635-53777704': ('chr5:174007939-174009432', 'chr5:130752381-130753830'), 'chr16:53925426-53926800': ('chr16:54315183-54322626', 'chr16:53925426-53926800'), 'chr16:53979176-53980851': ('chr16:54315183-54322626', 'chr16:53925426-53926800'), 'chr16:54315183-54322626': ('chr16:53088027-53090207', 'chr16:86527974-86528787')}, '/input_dir/hES_hsa_400kb_xyz.bed', 'http://genome.ucsc.edu/cgi-bin/hgTracks?hgS_doOt...Name=euplotidUCSC&hgS_otherUserSessionName=primed')
            kwargs = {}
        132 
        133     def __len__(self):
        134         return self._size
        135 
    
    ...........................................................................
    /root/Euplotid/src/buildAnnotateINs.py in build_IN_levels(target_node_name='chr16:53979176-53980851', dendogram_com=[{'chr10:100002276-100003364': 108, 'chr10:100021906-100023080': 108, 'chr10:100067964-100070284': 108, 'chr10:100091165-100092468': 108, 'chr10:100174343-100175470': 108, 'chr10:100183011-100184125': 1528, 'chr10:100185094-100186268': 1528, 'chr10:100205605-100207134': 7231, 'chr10:100226542-100228626': 7231, 'chr10:100423690-100425020': 1528, ...}, {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7, 8: 8, 9: 9, ...}, {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7, 8: 8, 9: 9, ...}], dna_int_graph=<networkx.classes.graph.Graph object>, genome_fa='/input_dir/hg19.fa', chain_file='/input_dir/hg19ToHg38.over.chain.gz', chain_file2='/input_dir/hg38ToHg19.over.chain.gz', out_dir='/output_dir//primed/example/', run_name='example', json_draft_name='/output_dir//primed/example/example_global_graph.json', TF_RBP_ids='/input_dir/TF_RBP_human.ids', open_peaks_file='open_peaks_nodes.bed', picked_tissue='Adipose_Visceral_Omentum', ctcf2bounds={'chr16:53088027-53090207': ('chr16:53088027-53090207', 'chr16:86527974-86528787'), 'chr16:53132274-53134474': ('chr16:53088027-53090207', 'chr16:86527974-86528787'), 'chr16:53193391-53194589': ('chr16:53088027-53090207', 'chr16:53193391-53194589'), 'chr16:53241256-53242869': ('chr16:53055849-53057134', 'chr16:52101824-52104715'), 'chr16:53254174-53255315': ('chr16:53088027-53090207', 'chr16:53254174-53255315'), 'chr16:53776635-53777704': ('chr5:174007939-174009432', 'chr5:130752381-130753830'), 'chr16:53925426-53926800': ('chr16:54315183-54322626', 'chr16:53925426-53926800'), 'chr16:53979176-53980851': ('chr16:54315183-54322626', 'chr16:53925426-53926800'), 'chr16:54315183-54322626': ('chr16:53088027-53090207', 'chr16:86527974-86528787')}, hsa_3D_loc='/input_dir/hES_hsa_400kb_xyz.bed', ucsc_session='http://genome.ucsc.edu/cgi-bin/hgTracks?hgS_doOt...Name=euplotidUCSC&hgS_otherUserSessionName=primed')
        776         add_snp = True
        777     else:
        778         add_snp = False
        779     source_community_graph = annotate_IN(target_node_name, dendogram_com, level, dna_int_graph, add_open, add_snp,
        780                                                   homo_gen, chain_file, chain_file2, TF_RBP_ids, open_peaks_file, picked_tissue, ctcf2bounds, comm_open_out,
    --> 781                                                   comm_var_out, comm_nodes_out, comm_edges_out, comm_edges_out_bed, ucsc_session_in,out_dir)
            comm_var_out = <_io.TextIOWrapper name='/output_dir//primed/exa...y_variants_lvl_1.bed' mode='a+' encoding='UTF-8'>
            comm_nodes_out = <_io.TextIOWrapper name='/output_dir//primed/exa...nity_nodes_lvl_1.bed' mode='a+' encoding='UTF-8'>
            comm_edges_out = <_io.TextIOWrapper name='/output_dir//primed/exa...ty_edges_lvl_1.washu' mode='a+' encoding='UTF-8'>
            comm_edges_out_bed = <_io.TextIOWrapper name='/output_dir//primed/exa...nity_edges_lvl_1.bed' mode='a+' encoding='UTF-8'>
            ucsc_session_in = 'http://genome.ucsc.edu/cgi-bin/hgTracks?hgS_doOt...ssionName=primed&position=chr16:54668601-53806737'
            out_dir = '/output_dir//primed/example/'
        782     #level 1 
        783     if len(dendogram_com) > 1:
        784         level = 1
        785         add_open = False
    
    ...........................................................................
    /root/Euplotid/src/buildAnnotateINs.py in annotate_IN(target_node_name='chr16:53979176-53980851', dendogram_com=[{'chr10:100002276-100003364': 108, 'chr10:100021906-100023080': 108, 'chr10:100067964-100070284': 108, 'chr10:100091165-100092468': 108, 'chr10:100174343-100175470': 108, 'chr10:100183011-100184125': 1528, 'chr10:100185094-100186268': 1528, 'chr10:100205605-100207134': 7231, 'chr10:100226542-100228626': 7231, 'chr10:100423690-100425020': 1528, ...}, {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7, 8: 8, 9: 9, ...}, {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7, 8: 8, 9: 9, ...}], level=0, dna_int_graph=<networkx.classes.graph.Graph object>, add_open=True, add_snp=True, homo_gen=<pysam.libcfaidx.FastaFile object>, chain_file='/input_dir/hg19ToHg38.over.chain.gz', chain_file2='/input_dir/hg38ToHg19.over.chain.gz', TF_RBP_ids='/input_dir/TF_RBP_human.ids', open_peaks_file='open_peaks_nodes.bed', picked_tissue='Adipose_Visceral_Omentum', ctcf2bounds={'chr16:53088027-53090207': ('chr16:53088027-53090207', 'chr16:86527974-86528787'), 'chr16:53132274-53134474': ('chr16:53088027-53090207', 'chr16:86527974-86528787'), 'chr16:53193391-53194589': ('chr16:53088027-53090207', 'chr16:53193391-53194589'), 'chr16:53241256-53242869': ('chr16:53055849-53057134', 'chr16:52101824-52104715'), 'chr16:53254174-53255315': ('chr16:53088027-53090207', 'chr16:53254174-53255315'), 'chr16:53776635-53777704': ('chr5:174007939-174009432', 'chr5:130752381-130753830'), 'chr16:53925426-53926800': ('chr16:54315183-54322626', 'chr16:53925426-53926800'), 'chr16:53979176-53980851': ('chr16:54315183-54322626', 'chr16:53925426-53926800'), 'chr16:54315183-54322626': ('chr16:53088027-53090207', 'chr16:86527974-86528787')}, comm_open_out=<_io.TextIOWrapper name='/output_dir//primed/exa...en_regions_lvl_1.bed' mode='a+' encoding='UTF-8'>, comm_var_out=<_io.TextIOWrapper name='/output_dir//primed/exa...y_variants_lvl_1.bed' mode='a+' encoding='UTF-8'>, comm_nodes_out=<_io.TextIOWrapper name='/output_dir//primed/exa...nity_nodes_lvl_1.bed' mode='a+' encoding='UTF-8'>, comm_edges_out=<_io.TextIOWrapper name='/output_dir//primed/exa...ty_edges_lvl_1.washu' mode='a+' encoding='UTF-8'>, comm_edges_out_bed=<_io.TextIOWrapper name='/output_dir//primed/exa...nity_edges_lvl_1.bed' mode='a+' encoding='UTF-8'>, ucsc_session_in='http://genome.ucsc.edu/cgi-bin/hgTracks?hgS_doOt...ssionName=primed&position=chr16:54668601-53806737', out_dir='/output_dir//primed/example/')
        668         print("Annotating Open Regions")
        669         source_community_graph = add_openRegions_predict(source_community_graph, homo_gen, chain_file, chain_file2, TF_RBP_ids, open_peaks_file, target_node_name, level, dict())
        670         sys.stdout.flush()
        671         if add_snp:
        672             print("Annotating SNPs/CNVs")
    --> 673             source_community_graph = add_variants_predict(source_community_graph, homo_gen, chain_file, TF_RBP_ids, picked_tissue, target_node_name,out_dir)
            source_community_graph = <networkx.classes.graph.Graph object>
            homo_gen = <pysam.libcfaidx.FastaFile object>
            chain_file = '/input_dir/hg19ToHg38.over.chain.gz'
            TF_RBP_ids = '/input_dir/TF_RBP_human.ids'
            picked_tissue = 'Adipose_Visceral_Omentum'
            target_node_name = 'chr16:53979176-53980851'
            out_dir = '/output_dir//primed/example/'
        674             sys.stdout.flush()
        675     
        676     #Color Insulated Neighborhood according to node attributes, make mid distance to target node
        677 
    
    ...........................................................................
    /root/Euplotid/src/buildAnnotateINs.py in add_variants_predict(G=<networkx.classes.graph.Graph object>, homo_gen=<pysam.libcfaidx.FastaFile object>, chain_file='/input_dir/hg19ToHg38.over.chain.gz', TF_RBP_ids='/input_dir/TF_RBP_human.ids', tissue_type='Adipose_Visceral_Omentum', target_node_name='chr16:53979176-53980851', out_dir='/output_dir//primed/example/')
        297     #predict effect of all variation added to graph using DeepBind and Basset
        298     if add_any_SNPs:
        299         print("Deepbind SNP prediction of: " + target_node_name)
        300         G = deepbind_predict_SNPs(G, homo_gen, TF_RBP_ids,target_node_name)
        301         print("Basset SNP prediction of: " + target_node_name)
    --> 302         G = add_basset_sad_sat(G, homo_gen, target_node_name,out_dir)
            G = <networkx.classes.graph.Graph object>
            homo_gen = <pysam.libcfaidx.FastaFile object>
            target_node_name = 'chr16:53979176-53980851'
            out_dir = '/output_dir//primed/example/'
        303     try:
        304         os.remove(basset_name_out + "open_variants.vcf")
        305     except OSError:
        306         pass
    
    ...........................................................................
    /root/Euplotid/src/buildAnnotateINs.py in add_basset_sad_sat(G=<networkx.classes.graph.Graph object>, homo_gen=<pysam.libcfaidx.FastaFile object>, target_node_name='chr16:53979176-53980851', out_dir='/output_dir//primed/example/')
        401             count += 1
        402     cmd = ("/root/Basset/src/basset_sad.py -f " + (homo_gen.filename).decode("utf-8") + " -l 600 -o " + basset_out_name + " -t " + targets_file + " " +  model_file + " " + basset_out_name + "/open_variants.vcf")
        403     print(cmd)
        404     os.system(cmd)
        405     #Pick SNP maximizing for sum of Delta SAD across all cell types in trained model
    --> 406     sad_table = pd.read_table(basset_out_name + "/sad_table.txt", delim_whitespace=True, header = 0)
            sad_table = undefined
            basset_out_name = '/output_dir//primed/example//chr16_53979176_53980851_basset'
        407     sad_table_dense = sad_table.pivot_table(index='rsid', columns='target', values='pred')
        408     serial_sad_table = pickle.dumps(sad_table_dense, protocol=0)
        409     G.graph["sad_table"] = serial_sad_table
        410     sad_table["sad_abs"] = abs(sad_table["pred"])
    
    ...........................................................................
    /opt/conda/lib/python3.5/site-packages/pandas/io/parsers.py in parser_f(filepath_or_buffer='/output_dir//primed/example//chr16_53979176_53980851_basset/sad_table.txt', sep='\t', delimiter='\t', header=0, names=None, index_col=None, usecols=None, squeeze=False, prefix=None, mangle_dupe_cols=True, dtype=None, engine='c', converters=None, true_values=None, false_values=None, skipinitialspace=False, skiprows=None, nrows=None, na_values=None, keep_default_na=True, na_filter=True, verbose=False, skip_blank_lines=True, parse_dates=False, infer_datetime_format=False, keep_date_col=False, date_parser=None, dayfirst=False, iterator=False, chunksize=None, compression='infer', thousands=None, decimal=b'.', lineterminator=None, quotechar='"', quoting=0, escapechar=None, comment=None, encoding=None, dialect=None, tupleize_cols=False, error_bad_lines=True, warn_bad_lines=True, skipfooter=0, skip_footer=0, doublequote=True, delim_whitespace=True, as_recarray=False, compact_ints=False, use_unsigned=False, low_memory=True, buffer_lines=None, memory_map=False, float_precision=None)
        650                     mangle_dupe_cols=mangle_dupe_cols,
        651                     tupleize_cols=tupleize_cols,
        652                     infer_datetime_format=infer_datetime_format,
        653                     skip_blank_lines=skip_blank_lines)
        654 
    --> 655         return _read(filepath_or_buffer, kwds)
            filepath_or_buffer = '/output_dir//primed/example//chr16_53979176_53980851_basset/sad_table.txt'
            kwds = {'as_recarray': False, 'buffer_lines': None, 'chunksize': None, 'comment': None, 'compact_ints': False, 'compression': None, 'converters': None, 'date_parser': None, 'dayfirst': False, 'decimal': b'.', ...}
        656 
        657     parser_f.__name__ = name
        658 
        659     return parser_f
    
    ...........................................................................
    /opt/conda/lib/python3.5/site-packages/pandas/io/parsers.py in _read(filepath_or_buffer='/output_dir//primed/example//chr16_53979176_53980851_basset/sad_table.txt', kwds={'as_recarray': False, 'buffer_lines': None, 'chunksize': None, 'comment': None, 'compact_ints': False, 'compression': None, 'converters': None, 'date_parser': None, 'dayfirst': False, 'decimal': b'.', ...})
        400     iterator = kwds.get('iterator', False)
        401     chunksize = _validate_integer('chunksize', kwds.get('chunksize', None), 1)
        402     nrows = _validate_integer('nrows', kwds.get('nrows', None))
        403 
        404     # Create the parser.
    --> 405     parser = TextFileReader(filepath_or_buffer, **kwds)
            parser = undefined
            filepath_or_buffer = '/output_dir//primed/example//chr16_53979176_53980851_basset/sad_table.txt'
            kwds = {'as_recarray': False, 'buffer_lines': None, 'chunksize': None, 'comment': None, 'compact_ints': False, 'compression': None, 'converters': None, 'date_parser': None, 'dayfirst': False, 'decimal': b'.', ...}
        406 
        407     if chunksize or iterator:
        408         return parser
        409 
    
    ...........................................................................
    /opt/conda/lib/python3.5/site-packages/pandas/io/parsers.py in __init__(self=<pandas.io.parsers.TextFileReader object>, f='/output_dir//primed/example//chr16_53979176_53980851_basset/sad_table.txt', engine='c', **kwds={'as_recarray': False, 'buffer_lines': None, 'chunksize': None, 'comment': None, 'compact_ints': False, 'compression': None, 'converters': None, 'date_parser': None, 'dayfirst': False, 'decimal': b'.', ...})
        757         # might mutate self.engine
        758         self.options, self.engine = self._clean_options(options, engine)
        759         if 'has_index_names' in kwds:
        760             self.options['has_index_names'] = kwds['has_index_names']
        761 
    --> 762         self._make_engine(self.engine)
            self._make_engine = <bound method TextFileReader._make_engine of <pandas.io.parsers.TextFileReader object>>
            self.engine = 'c'
        763 
        764     def close(self):
        765         self._engine.close()
        766 
    
    ...........................................................................
    /opt/conda/lib/python3.5/site-packages/pandas/io/parsers.py in _make_engine(self=<pandas.io.parsers.TextFileReader object>, engine='c')
        961             self.close()
        962             raise
        963 
        964     def _make_engine(self, engine='c'):
        965         if engine == 'c':
    --> 966             self._engine = CParserWrapper(self.f, **self.options)
            self._engine = None
            self.f = '/output_dir//primed/example//chr16_53979176_53980851_basset/sad_table.txt'
            self.options = {'as_recarray': False, 'buffer_lines': None, 'comment': None, 'compact_ints': False, 'compression': None, 'converters': {}, 'date_parser': None, 'dayfirst': False, 'decimal': b'.', 'delim_whitespace': True, ...}
        967         else:
        968             if engine == 'python':
        969                 klass = PythonParser
        970             elif engine == 'python-fwf':
    
    ...........................................................................
    /opt/conda/lib/python3.5/site-packages/pandas/io/parsers.py in __init__(self=<pandas.io.parsers.CParserWrapper object>, src='/output_dir//primed/example//chr16_53979176_53980851_basset/sad_table.txt', **kwds={'allow_leading_cols': True, 'as_recarray': False, 'buffer_lines': None, 'comment': None, 'compact_ints': False, 'compression': None, 'converters': {}, 'decimal': b'.', 'delim_whitespace': True, 'delimiter': '\t', ...})
       1577             kwds['encoding'] = 'utf-8'
       1578 
       1579         # #2442
       1580         kwds['allow_leading_cols'] = self.index_col is not False
       1581 
    -> 1582         self._reader = parsers.TextReader(src, **kwds)
            self._reader = undefined
            src = '/output_dir//primed/example//chr16_53979176_53980851_basset/sad_table.txt'
            kwds = {'allow_leading_cols': True, 'as_recarray': False, 'buffer_lines': None, 'comment': None, 'compact_ints': False, 'compression': None, 'converters': {}, 'decimal': b'.', 'delim_whitespace': True, 'delimiter': '\t', ...}
       1583 
       1584         # XXX
       1585         self.usecols, self.usecols_dtype = _validate_usecols_arg(
       1586             self._reader.usecols)
    
    ...........................................................................
    /opt/conda/lib/python3.5/site-packages/pandas/_libs/parsers.cpython-35m-x86_64-linux-gnu.so in pandas._libs.parsers.TextReader.__cinit__ (pandas/_libs/parsers.c:4209)()
    
    ...........................................................................
    /opt/conda/lib/python3.5/site-packages/pandas/_libs/parsers.cpython-35m-x86_64-linux-gnu.so in pandas._libs.parsers.TextReader._setup_parser_source (pandas/_libs/parsers.c:8873)()
    
    FileNotFoundError: File b'/output_dir//primed/example//chr16_53979176_53980851_basset/sad_table.txt' does not exist
    ___________________________________________________________________________
    
    During handling of the above exception, another exception occurred:
    
    Traceback (most recent call last):
      File "/root/Euplotid/src/buildAnnotateINs.py", line 1195, in <module>
        TF_RBP_ids, open_peaks_file, picked_tissue, ctcf2bounds, hsa_3D_loc, ucsc_session) for target_node_name in final_in_nodes)
      File "/opt/conda/lib/python3.5/site-packages/joblib/parallel.py", line 789, in __call__
        self.retrieve()
      File "/opt/conda/lib/python3.5/site-packages/joblib/parallel.py", line 740, in retrieve
        raise exception
    joblib.my_exceptions.JoblibFileNotFoundError: JoblibFileNotFoundError
    ___________________________________________________________________________
    Multiprocessing exception:
    ...........................................................................
    /root/Euplotid/src/buildAnnotateINs.py in <module>()
       1190     os.system("closest-features --closest --dist " + open_peaks + " all_nodes_sorted.bed > " + open_peaks_file)
       1191     
       1192     num_cores = multiprocessing.cpu_count()
       1193 
       1194     results = Parallel(n_jobs=num_cores,temp_folder=out_dir) (delayed(build_IN_levels) (target_node_name, dendogram_com, dna_int_graph, genome_fa, chain_file, chain_file2, out_dir, run_name, json_draft_name,
    -> 1195                                                                     TF_RBP_ids, open_peaks_file, picked_tissue, ctcf2bounds, hsa_3D_loc, ucsc_session) for target_node_name in final_in_nodes)
       1196 
       1197     #Add check if succeded?
       1198     #Annotate draft genome    
       1199     json_comm = open(json_draft_name,"r")
    
    ...........................................................................
    /opt/conda/lib/python3.5/site-packages/joblib/parallel.py in __call__(self=Parallel(n_jobs=6), iterable=<generator object <genexpr>>)
        784             if pre_dispatch == "all" or n_jobs == 1:
        785                 # The iterable was consumed all at once by the above for loop.
        786                 # No need to wait for async callbacks to trigger to
        787                 # consumption.
        788                 self._iterating = False
    --> 789             self.retrieve()
            self.retrieve = <bound method Parallel.retrieve of Parallel(n_jobs=6)>
        790             # Make sure that we get a last message telling us we are done
        791             elapsed_time = time.time() - self._start_time
        792             self._print('Done %3i out of %3i | elapsed: %s finished',
        793                         (len(self._output), len(self._output),
    
    ---------------------------------------------------------------------------
    Sub-process traceback:
    ---------------------------------------------------------------------------
    FileNotFoundError                                  Wed May 31 03:20:16 2017
    PID: 3714                               Python 3.5.3: /opt/conda/bin/python
    ...........................................................................
    /opt/conda/lib/python3.5/site-packages/joblib/parallel.py in __call__(self=<joblib.parallel.BatchedCalls object>)
        126     def __init__(self, iterator_slice):
        127         self.items = list(iterator_slice)
        128         self._size = len(self.items)
        129 
        130     def __call__(self):
    --> 131         return [func(*args, **kwargs) for func, args, kwargs in self.items]
            self.items = [(<function build_IN_levels>, ('chr16:53979176-53980851', [{'chr10:100002276-100003364': 108, 'chr10:100021906-100023080': 108, 'chr10:100067964-100070284': 108, 'chr10:100091165-100092468': 108, 'chr10:100174343-100175470': 108, 'chr10:100183011-100184125': 1475, 'chr10:100185094-100186268': 1475, 'chr10:100205605-100207134': 1475, 'chr10:100226542-100228626': 1475, 'chr10:100423690-100425020': 1475, ...}, {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7, 8: 8, 9: 9, ...}, {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7, 8: 8, 9: 9, ...}, {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7, 8: 8, 9: 9, ...}], <networkx.classes.graph.Graph object>, '/input_dir/hg19.fa', '/input_dir/hg19ToHg38.over.chain.gz', '/input_dir/hg38ToHg19.over.chain.gz', '/output_dir//primed/example/', 'example', '/output_dir//primed/example/example_global_graph.json', '/input_dir/TF_RBP_human.ids', 'open_peaks_nodes.bed', 'Adipose_Visceral_Omentum', {'chr16:53088027-53090207': ('chr16:53088027-53090207', 'chr16:86527974-86528787'), 'chr16:53132274-53134474': ('chr16:53088027-53090207', 'chr16:86527974-86528787'), 'chr16:53193391-53194589': ('chr16:53088027-53090207', 'chr16:53193391-53194589'), 'chr16:53241256-53242869': ('chr16:53055849-53057134', 'chr16:52101824-52104715'), 'chr16:53254174-53255315': ('chr16:53088027-53090207', 'chr16:53254174-53255315'), 'chr16:53776635-53777704': ('chr5:174007939-174009432', 'chr5:130752381-130753830'), 'chr16:53925426-53926800': ('chr16:54315183-54322626', 'chr16:53925426-53926800'), 'chr16:53979176-53980851': ('chr16:54315183-54322626', 'chr16:53925426-53926800'), 'chr16:54315183-54322626': ('chr16:53088027-53090207', 'chr16:86527974-86528787')}, '/input_dir/hES_hsa_400kb_xyz.bed', 'http://genome.ucsc.edu/cgi-bin/hgTracks?hgS_doOt...Name=euplotidUCSC&hgS_otherUserSessionName=primed'), {})]
        132 
        133     def __len__(self):
        134         return self._size
        135 
    
    ...........................................................................
    /opt/conda/lib/python3.5/site-packages/joblib/parallel.py in <listcomp>(.0=<list_iterator object>)
        126     def __init__(self, iterator_slice):
        127         self.items = list(iterator_slice)
        128         self._size = len(self.items)
        129 
        130     def __call__(self):
    --> 131         return [func(*args, **kwargs) for func, args, kwargs in self.items]
            func = <function build_IN_levels>
            args = ('chr16:53979176-53980851', [{'chr10:100002276-100003364': 108, 'chr10:100021906-100023080': 108, 'chr10:100067964-100070284': 108, 'chr10:100091165-100092468': 108, 'chr10:100174343-100175470': 108, 'chr10:100183011-100184125': 1475, 'chr10:100185094-100186268': 1475, 'chr10:100205605-100207134': 1475, 'chr10:100226542-100228626': 1475, 'chr10:100423690-100425020': 1475, ...}, {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7, 8: 8, 9: 9, ...}, {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7, 8: 8, 9: 9, ...}, {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7, 8: 8, 9: 9, ...}], <networkx.classes.graph.Graph object>, '/input_dir/hg19.fa', '/input_dir/hg19ToHg38.over.chain.gz', '/input_dir/hg38ToHg19.over.chain.gz', '/output_dir//primed/example/', 'example', '/output_dir//primed/example/example_global_graph.json', '/input_dir/TF_RBP_human.ids', 'open_peaks_nodes.bed', 'Adipose_Visceral_Omentum', {'chr16:53088027-53090207': ('chr16:53088027-53090207', 'chr16:86527974-86528787'), 'chr16:53132274-53134474': ('chr16:53088027-53090207', 'chr16:86527974-86528787'), 'chr16:53193391-53194589': ('chr16:53088027-53090207', 'chr16:53193391-53194589'), 'chr16:53241256-53242869': ('chr16:53055849-53057134', 'chr16:52101824-52104715'), 'chr16:53254174-53255315': ('chr16:53088027-53090207', 'chr16:53254174-53255315'), 'chr16:53776635-53777704': ('chr5:174007939-174009432', 'chr5:130752381-130753830'), 'chr16:53925426-53926800': ('chr16:54315183-54322626', 'chr16:53925426-53926800'), 'chr16:53979176-53980851': ('chr16:54315183-54322626', 'chr16:53925426-53926800'), 'chr16:54315183-54322626': ('chr16:53088027-53090207', 'chr16:86527974-86528787')}, '/input_dir/hES_hsa_400kb_xyz.bed', 'http://genome.ucsc.edu/cgi-bin/hgTracks?hgS_doOt...Name=euplotidUCSC&hgS_otherUserSessionName=primed')
            kwargs = {}
        132 
        133     def __len__(self):
        134         return self._size
        135 
    
    ...........................................................................
    /root/Euplotid/src/buildAnnotateINs.py in build_IN_levels(target_node_name='chr16:53979176-53980851', dendogram_com=[{'chr10:100002276-100003364': 108, 'chr10:100021906-100023080': 108, 'chr10:100067964-100070284': 108, 'chr10:100091165-100092468': 108, 'chr10:100174343-100175470': 108, 'chr10:100183011-100184125': 1528, 'chr10:100185094-100186268': 1528, 'chr10:100205605-100207134': 7231, 'chr10:100226542-100228626': 7231, 'chr10:100423690-100425020': 1528, ...}, {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7, 8: 8, 9: 9, ...}, {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7, 8: 8, 9: 9, ...}], dna_int_graph=<networkx.classes.graph.Graph object>, genome_fa='/input_dir/hg19.fa', chain_file='/input_dir/hg19ToHg38.over.chain.gz', chain_file2='/input_dir/hg38ToHg19.over.chain.gz', out_dir='/output_dir//primed/example/', run_name='example', json_draft_name='/output_dir//primed/example/example_global_graph.json', TF_RBP_ids='/input_dir/TF_RBP_human.ids', open_peaks_file='open_peaks_nodes.bed', picked_tissue='Adipose_Visceral_Omentum', ctcf2bounds={'chr16:53088027-53090207': ('chr16:53088027-53090207', 'chr16:86527974-86528787'), 'chr16:53132274-53134474': ('chr16:53088027-53090207', 'chr16:86527974-86528787'), 'chr16:53193391-53194589': ('chr16:53088027-53090207', 'chr16:53193391-53194589'), 'chr16:53241256-53242869': ('chr16:53055849-53057134', 'chr16:52101824-52104715'), 'chr16:53254174-53255315': ('chr16:53088027-53090207', 'chr16:53254174-53255315'), 'chr16:53776635-53777704': ('chr5:174007939-174009432', 'chr5:130752381-130753830'), 'chr16:53925426-53926800': ('chr16:54315183-54322626', 'chr16:53925426-53926800'), 'chr16:53979176-53980851': ('chr16:54315183-54322626', 'chr16:53925426-53926800'), 'chr16:54315183-54322626': ('chr16:53088027-53090207', 'chr16:86527974-86528787')}, hsa_3D_loc='/input_dir/hES_hsa_400kb_xyz.bed', ucsc_session='http://genome.ucsc.edu/cgi-bin/hgTracks?hgS_doOt...Name=euplotidUCSC&hgS_otherUserSessionName=primed')
        776         add_snp = True
        777     else:
        778         add_snp = False
        779     source_community_graph = annotate_IN(target_node_name, dendogram_com, level, dna_int_graph, add_open, add_snp,
        780                                                   homo_gen, chain_file, chain_file2, TF_RBP_ids, open_peaks_file, picked_tissue, ctcf2bounds, comm_open_out,
    --> 781                                                   comm_var_out, comm_nodes_out, comm_edges_out, comm_edges_out_bed, ucsc_session_in,out_dir)
            comm_var_out = <_io.TextIOWrapper name='/output_dir//primed/exa...y_variants_lvl_1.bed' mode='a+' encoding='UTF-8'>
            comm_nodes_out = <_io.TextIOWrapper name='/output_dir//primed/exa...nity_nodes_lvl_1.bed' mode='a+' encoding='UTF-8'>
            comm_edges_out = <_io.TextIOWrapper name='/output_dir//primed/exa...ty_edges_lvl_1.washu' mode='a+' encoding='UTF-8'>
            comm_edges_out_bed = <_io.TextIOWrapper name='/output_dir//primed/exa...nity_edges_lvl_1.bed' mode='a+' encoding='UTF-8'>
            ucsc_session_in = 'http://genome.ucsc.edu/cgi-bin/hgTracks?hgS_doOt...ssionName=primed&position=chr16:54668601-53806737'
            out_dir = '/output_dir//primed/example/'
        782     #level 1 
        783     if len(dendogram_com) > 1:
        784         level = 1
        785         add_open = False
    
    ...........................................................................
    /root/Euplotid/src/buildAnnotateINs.py in annotate_IN(target_node_name='chr16:53979176-53980851', dendogram_com=[{'chr10:100002276-100003364': 108, 'chr10:100021906-100023080': 108, 'chr10:100067964-100070284': 108, 'chr10:100091165-100092468': 108, 'chr10:100174343-100175470': 108, 'chr10:100183011-100184125': 1528, 'chr10:100185094-100186268': 1528, 'chr10:100205605-100207134': 7231, 'chr10:100226542-100228626': 7231, 'chr10:100423690-100425020': 1528, ...}, {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7, 8: 8, 9: 9, ...}, {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7, 8: 8, 9: 9, ...}], level=0, dna_int_graph=<networkx.classes.graph.Graph object>, add_open=True, add_snp=True, homo_gen=<pysam.libcfaidx.FastaFile object>, chain_file='/input_dir/hg19ToHg38.over.chain.gz', chain_file2='/input_dir/hg38ToHg19.over.chain.gz', TF_RBP_ids='/input_dir/TF_RBP_human.ids', open_peaks_file='open_peaks_nodes.bed', picked_tissue='Adipose_Visceral_Omentum', ctcf2bounds={'chr16:53088027-53090207': ('chr16:53088027-53090207', 'chr16:86527974-86528787'), 'chr16:53132274-53134474': ('chr16:53088027-53090207', 'chr16:86527974-86528787'), 'chr16:53193391-53194589': ('chr16:53088027-53090207', 'chr16:53193391-53194589'), 'chr16:53241256-53242869': ('chr16:53055849-53057134', 'chr16:52101824-52104715'), 'chr16:53254174-53255315': ('chr16:53088027-53090207', 'chr16:53254174-53255315'), 'chr16:53776635-53777704': ('chr5:174007939-174009432', 'chr5:130752381-130753830'), 'chr16:53925426-53926800': ('chr16:54315183-54322626', 'chr16:53925426-53926800'), 'chr16:53979176-53980851': ('chr16:54315183-54322626', 'chr16:53925426-53926800'), 'chr16:54315183-54322626': ('chr16:53088027-53090207', 'chr16:86527974-86528787')}, comm_open_out=<_io.TextIOWrapper name='/output_dir//primed/exa...en_regions_lvl_1.bed' mode='a+' encoding='UTF-8'>, comm_var_out=<_io.TextIOWrapper name='/output_dir//primed/exa...y_variants_lvl_1.bed' mode='a+' encoding='UTF-8'>, comm_nodes_out=<_io.TextIOWrapper name='/output_dir//primed/exa...nity_nodes_lvl_1.bed' mode='a+' encoding='UTF-8'>, comm_edges_out=<_io.TextIOWrapper name='/output_dir//primed/exa...ty_edges_lvl_1.washu' mode='a+' encoding='UTF-8'>, comm_edges_out_bed=<_io.TextIOWrapper name='/output_dir//primed/exa...nity_edges_lvl_1.bed' mode='a+' encoding='UTF-8'>, ucsc_session_in='http://genome.ucsc.edu/cgi-bin/hgTracks?hgS_doOt...ssionName=primed&position=chr16:54668601-53806737', out_dir='/output_dir//primed/example/')
        668         print("Annotating Open Regions")
        669         source_community_graph = add_openRegions_predict(source_community_graph, homo_gen, chain_file, chain_file2, TF_RBP_ids, open_peaks_file, target_node_name, level, dict())
        670         sys.stdout.flush()
        671         if add_snp:
        672             print("Annotating SNPs/CNVs")
    --> 673             source_community_graph = add_variants_predict(source_community_graph, homo_gen, chain_file, TF_RBP_ids, picked_tissue, target_node_name,out_dir)
            source_community_graph = <networkx.classes.graph.Graph object>
            homo_gen = <pysam.libcfaidx.FastaFile object>
            chain_file = '/input_dir/hg19ToHg38.over.chain.gz'
            TF_RBP_ids = '/input_dir/TF_RBP_human.ids'
            picked_tissue = 'Adipose_Visceral_Omentum'
            target_node_name = 'chr16:53979176-53980851'
            out_dir = '/output_dir//primed/example/'
        674             sys.stdout.flush()
        675     
        676     #Color Insulated Neighborhood according to node attributes, make mid distance to target node
        677 
    
    ...........................................................................
    /root/Euplotid/src/buildAnnotateINs.py in add_variants_predict(G=<networkx.classes.graph.Graph object>, homo_gen=<pysam.libcfaidx.FastaFile object>, chain_file='/input_dir/hg19ToHg38.over.chain.gz', TF_RBP_ids='/input_dir/TF_RBP_human.ids', tissue_type='Adipose_Visceral_Omentum', target_node_name='chr16:53979176-53980851', out_dir='/output_dir//primed/example/')
        297     #predict effect of all variation added to graph using DeepBind and Basset
        298     if add_any_SNPs:
        299         print("Deepbind SNP prediction of: " + target_node_name)
        300         G = deepbind_predict_SNPs(G, homo_gen, TF_RBP_ids,target_node_name)
        301         print("Basset SNP prediction of: " + target_node_name)
    --> 302         G = add_basset_sad_sat(G, homo_gen, target_node_name,out_dir)
            G = <networkx.classes.graph.Graph object>
            homo_gen = <pysam.libcfaidx.FastaFile object>
            target_node_name = 'chr16:53979176-53980851'
            out_dir = '/output_dir//primed/example/'
        303     try:
        304         os.remove(basset_name_out + "open_variants.vcf")
        305     except OSError:
        306         pass
    
    ...........................................................................
    /root/Euplotid/src/buildAnnotateINs.py in add_basset_sad_sat(G=<networkx.classes.graph.Graph object>, homo_gen=<pysam.libcfaidx.FastaFile object>, target_node_name='chr16:53979176-53980851', out_dir='/output_dir//primed/example/')
        401             count += 1
        402     cmd = ("/root/Basset/src/basset_sad.py -f " + (homo_gen.filename).decode("utf-8") + " -l 600 -o " + basset_out_name + " -t " + targets_file + " " +  model_file + " " + basset_out_name + "/open_variants.vcf")
        403     print(cmd)
        404     os.system(cmd)
        405     #Pick SNP maximizing for sum of Delta SAD across all cell types in trained model
    --> 406     sad_table = pd.read_table(basset_out_name + "/sad_table.txt", delim_whitespace=True, header = 0)
            sad_table = undefined
            basset_out_name = '/output_dir//primed/example//chr16_53979176_53980851_basset'
        407     sad_table_dense = sad_table.pivot_table(index='rsid', columns='target', values='pred')
        408     serial_sad_table = pickle.dumps(sad_table_dense, protocol=0)
        409     G.graph["sad_table"] = serial_sad_table
        410     sad_table["sad_abs"] = abs(sad_table["pred"])
    
    ...........................................................................
    /opt/conda/lib/python3.5/site-packages/pandas/io/parsers.py in parser_f(filepath_or_buffer='/output_dir//primed/example//chr16_53979176_53980851_basset/sad_table.txt', sep='\t', delimiter='\t', header=0, names=None, index_col=None, usecols=None, squeeze=False, prefix=None, mangle_dupe_cols=True, dtype=None, engine='c', converters=None, true_values=None, false_values=None, skipinitialspace=False, skiprows=None, nrows=None, na_values=None, keep_default_na=True, na_filter=True, verbose=False, skip_blank_lines=True, parse_dates=False, infer_datetime_format=False, keep_date_col=False, date_parser=None, dayfirst=False, iterator=False, chunksize=None, compression='infer', thousands=None, decimal=b'.', lineterminator=None, quotechar='"', quoting=0, escapechar=None, comment=None, encoding=None, dialect=None, tupleize_cols=False, error_bad_lines=True, warn_bad_lines=True, skipfooter=0, skip_footer=0, doublequote=True, delim_whitespace=True, as_recarray=False, compact_ints=False, use_unsigned=False, low_memory=True, buffer_lines=None, memory_map=False, float_precision=None)
        650                     mangle_dupe_cols=mangle_dupe_cols,
        651                     tupleize_cols=tupleize_cols,
        652                     infer_datetime_format=infer_datetime_format,
        653                     skip_blank_lines=skip_blank_lines)
        654 
    --> 655         return _read(filepath_or_buffer, kwds)
            filepath_or_buffer = '/output_dir//primed/example//chr16_53979176_53980851_basset/sad_table.txt'
            kwds = {'as_recarray': False, 'buffer_lines': None, 'chunksize': None, 'comment': None, 'compact_ints': False, 'compression': None, 'converters': None, 'date_parser': None, 'dayfirst': False, 'decimal': b'.', ...}
        656 
        657     parser_f.__name__ = name
        658 
        659     return parser_f
    
    ...........................................................................
    /opt/conda/lib/python3.5/site-packages/pandas/io/parsers.py in _read(filepath_or_buffer='/output_dir//primed/example//chr16_53979176_53980851_basset/sad_table.txt', kwds={'as_recarray': False, 'buffer_lines': None, 'chunksize': None, 'comment': None, 'compact_ints': False, 'compression': None, 'converters': None, 'date_parser': None, 'dayfirst': False, 'decimal': b'.', ...})
        400     iterator = kwds.get('iterator', False)
        401     chunksize = _validate_integer('chunksize', kwds.get('chunksize', None), 1)
        402     nrows = _validate_integer('nrows', kwds.get('nrows', None))
        403 
        404     # Create the parser.
    --> 405     parser = TextFileReader(filepath_or_buffer, **kwds)
            parser = undefined
            filepath_or_buffer = '/output_dir//primed/example//chr16_53979176_53980851_basset/sad_table.txt'
            kwds = {'as_recarray': False, 'buffer_lines': None, 'chunksize': None, 'comment': None, 'compact_ints': False, 'compression': None, 'converters': None, 'date_parser': None, 'dayfirst': False, 'decimal': b'.', ...}
        406 
        407     if chunksize or iterator:
        408         return parser
        409 
    
    ...........................................................................
    /opt/conda/lib/python3.5/site-packages/pandas/io/parsers.py in __init__(self=<pandas.io.parsers.TextFileReader object>, f='/output_dir//primed/example//chr16_53979176_53980851_basset/sad_table.txt', engine='c', **kwds={'as_recarray': False, 'buffer_lines': None, 'chunksize': None, 'comment': None, 'compact_ints': False, 'compression': None, 'converters': None, 'date_parser': None, 'dayfirst': False, 'decimal': b'.', ...})
        757         # might mutate self.engine
        758         self.options, self.engine = self._clean_options(options, engine)
        759         if 'has_index_names' in kwds:
        760             self.options['has_index_names'] = kwds['has_index_names']
        761 
    --> 762         self._make_engine(self.engine)
            self._make_engine = <bound method TextFileReader._make_engine of <pandas.io.parsers.TextFileReader object>>
            self.engine = 'c'
        763 
        764     def close(self):
        765         self._engine.close()
        766 
    
    ...........................................................................
    /opt/conda/lib/python3.5/site-packages/pandas/io/parsers.py in _make_engine(self=<pandas.io.parsers.TextFileReader object>, engine='c')
        961             self.close()
        962             raise
        963 
        964     def _make_engine(self, engine='c'):
        965         if engine == 'c':
    --> 966             self._engine = CParserWrapper(self.f, **self.options)
            self._engine = None
            self.f = '/output_dir//primed/example//chr16_53979176_53980851_basset/sad_table.txt'
            self.options = {'as_recarray': False, 'buffer_lines': None, 'comment': None, 'compact_ints': False, 'compression': None, 'converters': {}, 'date_parser': None, 'dayfirst': False, 'decimal': b'.', 'delim_whitespace': True, ...}
        967         else:
        968             if engine == 'python':
        969                 klass = PythonParser
        970             elif engine == 'python-fwf':
    
    ...........................................................................
    /opt/conda/lib/python3.5/site-packages/pandas/io/parsers.py in __init__(self=<pandas.io.parsers.CParserWrapper object>, src='/output_dir//primed/example//chr16_53979176_53980851_basset/sad_table.txt', **kwds={'allow_leading_cols': True, 'as_recarray': False, 'buffer_lines': None, 'comment': None, 'compact_ints': False, 'compression': None, 'converters': {}, 'decimal': b'.', 'delim_whitespace': True, 'delimiter': '\t', ...})
       1577             kwds['encoding'] = 'utf-8'
       1578 
       1579         # #2442
       1580         kwds['allow_leading_cols'] = self.index_col is not False
       1581 
    -> 1582         self._reader = parsers.TextReader(src, **kwds)
            self._reader = undefined
            src = '/output_dir//primed/example//chr16_53979176_53980851_basset/sad_table.txt'
            kwds = {'allow_leading_cols': True, 'as_recarray': False, 'buffer_lines': None, 'comment': None, 'compact_ints': False, 'compression': None, 'converters': {}, 'decimal': b'.', 'delim_whitespace': True, 'delimiter': '\t', ...}
       1583 
       1584         # XXX
       1585         self.usecols, self.usecols_dtype = _validate_usecols_arg(
       1586             self._reader.usecols)
    
    ...........................................................................
    /opt/conda/lib/python3.5/site-packages/pandas/_libs/parsers.cpython-35m-x86_64-linux-gnu.so in pandas._libs.parsers.TextReader.__cinit__ (pandas/_libs/parsers.c:4209)()
    
    ...........................................................................
    /opt/conda/lib/python3.5/site-packages/pandas/_libs/parsers.cpython-35m-x86_64-linux-gnu.so in pandas._libs.parsers.TextReader._setup_parser_source (pandas/_libs/parsers.c:8873)()
    
    FileNotFoundError: File b'/output_dir//primed/example//chr16_53979176_53980851_basset/sad_table.txt' does not exist
    ___________________________________________________________________________



```python
#Default Insulated Neighborhoods to 
!python generate_draft_JSONs_all_funcs.py -i /input_data/ -o /root/Euplotid/pipelines/example_output/output_IN_graphs/ \
    -c primed -q Adipose_Visceral_Omentum -n default 
!python generate_draft_JSONs_all_funcs.py -i /input_data/ -o /root/Euplotid/pipelines/example_output/output_IN_graphs/ \
    -c jurkatt -q Cells_EBV-transformed_lymphocytes -n default 
!python generate_draft_JSONs_all_funcs.py -i /input_data/ -o /root/Euplotid/pipelines/example_output/output_IN_graphs/ \
    -c neuron -q Brain_Cortex -n default 
!python generate_draft_JSONs_all_funcs.py -i /input_data/ -o /root/Euplotid/pipelines/example_output/output_IN_graphs/ \
    -c mesc -q Brain_Cortex -n default 
```

    /Users/dborgesr/anaconda/lib/python2.7/site-packages/matplotlib/font_manager.py:273: UserWarning:
    
    Matplotlib is building the font cache using fc-list. This may take a moment.
    

