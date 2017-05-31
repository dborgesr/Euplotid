
## Pipeline to find interconnected regions
This pipeline takes in an <a href="https://github.com/danielsday/origami/"> Origami </a> processed DNA-DNA interactions generated using [fq2ChIAInts](fq2ChIAInts.ipynb) and finds interconnected networks of nodes (communities) using default or picked resolution. Plotting functions are also included to visualize output  
Range for resolution is from [0,1] with 1e-4 being something close to Insulated Neighborhoods, higher resolution makes bigger communities. You can also pick the "level", we start at level 0 then merging level 0 communities we get level 1 communities, so on and so forth.


```python
# take a peek to see what DNA-DNA interactions are available
!ls /input_dir/*.bedpe -lt
#!conda install -y -c conda-forge  lua=5.3.3 
```

    -rw-r--r--    1 root     root        213315 May 22 02:32 [0;0m/input_dir/primed_.7_origami_chr6.bedpe[0m
    -rw-r--r--    1 root     root       4296019 May 21 01:36 [0;0m/input_dir/primed_.7_origami.bedpe[0m
    -rw-r--r--    1 root     root     270640140 Nov 27  2016 [0;0m/input_dir/Primed-SMC1-petcount.bedpe[0m
    -rw-r--r--    1 root     root        643953 Nov 27  2016 [0;0m/input_dir/NPC-SMC1-petcount-filtered.bedpe[0m
    -rw-r--r--    1 root     root       1798538 Nov 27  2016 [0;0m/input_dir/NPC_rep2_hg19.bedpe[0m
    -rw-r--r--    1 root     root       4501002 Nov 27  2016 [0;0m/input_dir/primed_n2.bedpe[0m
    -rw-r--r--    1 root     root       1288717 Nov 27  2016 [0;0m/input_dir/Neuron-SMC1-petcount-filtered.bedpe[0m
    -rw-r--r--    1 root     root       1357544 Nov 27  2016 [0;0m/input_dir/Neuron_SMC1_n2.bedpe[0m
    -rw-r--r--    1 root     root       4525609 Nov 27  2016 [0;0m/input_dir/Neuron_rep2_hg19.bedpe[0m
    -rw-r--r--    1 root     root       2346376 Nov 27  2016 [0;0m/input_dir/Naive-SMC1-petcount-filtered.bedpe[0m
    -rw-r--r--    1 root     root       9009631 Nov 27  2016 [0;0m/input_dir/mES-SMC1-merged-petcount-filtered.bedpe[0m
    -rw-r--r--    1 root     root       3317002 Nov 27  2016 [0;0m/input_dir/naive_n2.bedpe[0m
    -rw-r--r--    1 root     root       1188848 Nov 27  2016 [0;0m/input_dir/Jurkat-SMC1-petcount-filtered.bedpe[0m
    -rw-r--r--    1 root     root        886552 Nov 27  2016 [0;0m/input_dir/Jurkat_rep1and2_hg19_n3FDR.01.bedpe[0m
    -rw-r--r--    1 root     root       2055692 Nov 27  2016 [0;0m/input_dir/K562-RAD21-petcount-filtered.bedpe[0m
    -rw-r--r--    1 root     root        503605 Nov 27  2016 [0;0m/input_dir/mESC_SMC_HiChip.bedpe[0m
    -rw-r--r--    1 root     root        485934 Nov 27  2016 [0;0m/input_dir/GSE63525_GM12878_primary_replicate_HiCCUPS_rawcount.bedpe[0m
    -rw-r--r--    1 root     root        310031 Nov 27  2016 [0;0m/input_dir/GSE63525_K562_HiCCUPS_rawcount.bedpe[0m
    -rw-r--r--    1 root     root      11126437 Nov 27  2016 [0;0m/input_dir/GM12878-CTCF-filtered.bedpe[0m
    -rw-r--r--    1 root     root        521904 Nov 27  2016 [0;0m/input_dir/GM12878_SMC_HiChip.bedpe[0m
    -rw-r--r--    1 root     root        875350 Nov 27  2016 [0;0m/input_dir/GM12878_merged_.1_filt_mango_hg19.bedpe[0m
    -rw-r--r--    1 root     root        736521 Nov 27  2016 [0;0m/input_dir/ENCFF002EMO_GM12878_rep1_Rad21_Chia_hg19.bedpe[0m
    -rw-r--r--    1 root     root       1194896 Nov 27  2016 [0;0m/input_dir/Dowen_mesc_SMC1_chia_ints.bedpe[0m



```python
#Input DNA-DNA interactions and where should output text file of communities go?
dna_ints = "/input_dir/Dowen_mesc_SMC1_chia_ints.bedpe"
out_comms = "/output_dir/mES-SMC1_chiapet_commmunities.txt"
resolution = 1
pick_level = 0
```


```python
#import needed python packages
import networkx as nx, community, re, numpy as np, pandas as pd
from collections import Counter
import plotly
from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot
init_notebook_mode(connected=True)
import plotly.graph_objs as go
#Function to plot communities
def plot_comm(G):
    axis=dict(showbackground=False,
              showline=False,
              zeroline=False,
              showgrid=False,
              showticklabels=False,
              title="")
    layout = go.Layout(
             title=G.graph["bounds"],
             width=600,
             height=600,
             showlegend=False,
             xaxis=go.XAxis(axis),
             yaxis=go.YAxis(axis),
             scene=go.Scene(
             xaxis=go.XAxis(axis),
             yaxis=go.YAxis(axis),
             zaxis=go.ZAxis(axis)),
        margin=go.Margin(l=0,
                         r=0,
                         b=0,
                         t=50),
        hovermode="closest")
    all_position = nx.spring_layout(G,dim=2)  
    #make node trace
    traceN = go.Scatter(x=[], y=[], mode="markers", text=[],marker=go.Marker(color=[],size=[],opacity=[]))
    traceN["name"] = ""
    traceN["hoverinfo"] = "text"
    for all_node in G.nodes(data=True):
        text_node = all_node[0]
        traceN["text"].append(text_node.replace("\n","<br>"))
        traceN["marker"]["color"].append("grey")
        traceN["marker"]["size"].append(20)
        traceN["marker"]["opacity"].append(1)
        traceN["x"].append(all_position[all_node[0]][0])
        traceN["y"].append(all_position[all_node[0]][1])

    traceE = go.Scatter(x=[], y=[], mode="lines", hoverinfo = "none")
    traceE["name"] = ""
    traceE["line"]["width"] = 1
    for edge in G.edges(data=True):
        traceE["x"] += [all_position[edge[0]][0],all_position[edge[1]][0], None]
        traceE["y"] += [all_position[edge[0]][1],all_position[edge[1]][1], None]

    data = go.Data([traceE, traceN])
    fig = go.Figure(data=data, layout=layout)
    plotly.offline.iplot(fig)
```


<script>requirejs.config({paths: { 'plotly': ['https://cdn.plot.ly/plotly-latest.min']},});if(!window.Plotly) {{require(['plotly'],function(plotly) {window.Plotly=plotly;});}}</script>


## Louvain community finding
This chunk runs a <a href="https://en.wikipedia.org/wiki/Louvain_Modularity"> Louvain </a>  algorithm to partition the DNA-DNA interaction network into Communities. Output text file has the following format:
* First column is the start and end of the max / min of the community (linearly)
* Second column is the number of nodes that make up that community
* Third column and more are the member nodes


```python
#Make a graph to store DNA-DNA interactions
dna_int_graph = nx.Graph(style="filled")
#Load DNA interactions
with open(dna_ints) as dna_ints_iter:
    for dna_int in dna_ints_iter:
        arr = dna_int.split()
        x = str(arr[0]) + ":" + str(arr[1]) + "-" + str(arr[2])
        mid_x = (int(arr[1])+int(arr[2]))/2.0
        y = str(arr[3]) + ":" + str(arr[4]) + "-" + str(arr[5])
        mid_y = (int(arr[4])+int(arr[5]))/2.0
        if (arr[0] == arr[3]):
            dna_int_graph.add_edge(x,y,label=1, capacity = 1, weight=float(arr[6]))

#Run louvain algorithm at defined resolution
dendogram_com = community.generate_dendrogram(dna_int_graph, resolution = resolution)
#Pull the first level of partitions
dna_int_comm = community.partition_at_level(dendogram_com, pick_level)
size2comm = Counter(dna_int_comm.values())
size_comms = [size for size in size2comm.values()]
boundary2nodes = dict()
bound_list_num_nodes = list()
bound_list_names = list()
comm_num_nodes = list()
#iterate over dna_int_comm keys and add each start/end to set stored in dict of comm --> min and comm-->max
inv_map = {}
for k, v in dna_int_comm.items():
    inv_map[v] = inv_map.get(v, [])
    inv_map[v].append(k)
for comm, nodes in inv_map.items():
    #get max and min of neighborhoods
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
    comm_name = ref_node_arr[0] + ":" + str(min_comm) + "-" + str(max_comm+1)
    boundary2nodes[comm_name] = nodes
    bound_list_num_nodes.append(len(nodes))
    bound_list_names.append(comm_name)
#Output text file with first column as the boundary of community, second # nodes, third all nodes w/ tabs
with open(out_comms, "w") as out_txt_comm: 
    for key, value in boundary2nodes.items():
        out_txt_comm.write(key+"\t" + str(len(value)) + "\t")
        for node in value:
            out_txt_comm.write(node+"\t")
        out_txt_comm.write("\n")
```

# Find all Communities encompassing a DNA position
This allows you to pick a community to view as defined by a DNA position (ex. chr5:142275000), it will find ALL communities which encompass that location, print them, and plot the last one.


```python
#Where do you want to look?
dna_position = "chr1:13117222"
```


```python
#overlap with boundaries of defined communities and if it overlaps plot it
pick_split = re.split(r"[-:]",dna_position)
for key, value in boundary2nodes.items():
    bound_split = re.split(r"[-:]",key)
    if (pick_split[0] == bound_split[0]) and (int(bound_split[1])<int(pick_split[1])) and (int(bound_split[2])>int(pick_split[1])): 
        print("target: " + dna_position + "\tCommunity Overlapping: " 
              + key + "\tlinear_length: " + 
              str((int(bound_split[2])-int(bound_split[1]))/1000) + "Kb" +
              "\tnum_nodes: " + str(len(value)))
        pick_comm_graph = dna_int_graph.subgraph(value)
        pick_comm_graph.graph["bounds"] = key
plot_comm(pick_comm_graph)
```

    target: chr1:13117222	Community Overlapping: chr1:12854095-13129135	linear_length: 275.04Kb	num_nodes: 11



<div id="a275982b-6ec6-47d5-824a-c2ea5e574ff8" style="height: 600px; width: 600px;" class="plotly-graph-div"></div><script type="text/javascript">require(["plotly"], function(Plotly) { window.PLOTLYENV=window.PLOTLYENV || {};window.PLOTLYENV.BASE_URL="https://plot.ly";Plotly.newPlot("a275982b-6ec6-47d5-824a-c2ea5e574ff8", [{"mode": "lines", "hoverinfo": "none", "line": {"width": 1}, "x": [0.8457330444102182, 0.6472517033065207, null, 0.8457330444102182, 1.0, null, 0.88757081373272, 0.6472517033065207, null, 0.27609093398510565, 0.0, null, 0.27609093398510565, 0.4736678868127097, null, 0.27609093398510565, 0.0605507178772557, null, 0.3865383155915301, 0.4736678868127097, null, 0.6472517033065207, 0.4736678868127097, null, 0.6472517033065207, 0.8285873271603635, null, 0.4119271332613275, 0.4736678868127097, null], "name": "", "type": "scatter", "y": [0.5344176153380827, 0.3154454908125696, null, 0.5344176153380827, 0.7620378373639478, null, 0.2521932999515526, 0.3154454908125696, null, 0.342741894757011, 0.534157567282243, null, 0.342741894757011, 0.3267624048785176, null, 0.342741894757011, 0.2542842291367912, null, 0.021956346699537888, 0.3267624048785176, null, 0.3154454908125696, 0.3267624048785176, null, 0.3154454908125696, 0.0, null, 0.6320506156011613, 0.3267624048785176, null]}, {"mode": "markers", "hoverinfo": "text", "text": ["chr1:13121498-13123090", "chr1:13115518-13119048", "chr1:13093130-13095652", "chr1:13039584-13041392", "chr1:13049045-13050571", "chr1:13127929-13129134", "chr1:13041583-13043936", "chr1:12979738-12983484", "chr1:12854095-12856319", "chr1:13123202-13125339", "chr1:13107646-13108403"], "x": [0.8457330444102182, 1.0, 0.88757081373272, 0.27609093398510565, 0.3865383155915301, 0.0, 0.6472517033065207, 0.4119271332613275, 0.0605507178772557, 0.4736678868127097, 0.8285873271603635], "name": "", "type": "scatter", "marker": {"size": [20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20], "opacity": [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1], "color": ["grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey"]}, "y": [0.5344176153380827, 0.7620378373639478, 0.2521932999515526, 0.342741894757011, 0.021956346699537888, 0.534157567282243, 0.3154454908125696, 0.6320506156011613, 0.2542842291367912, 0.3267624048785176, 0.0]}], {"title": "chr1:12854095-13129135", "width": 600, "yaxis": {"title": "", "showgrid": false, "zeroline": false, "showbackground": false, "showticklabels": false, "showline": false}, "xaxis": {"title": "", "showgrid": false, "zeroline": false, "showbackground": false, "showticklabels": false, "showline": false}, "showlegend": false, "margin": {"r": 0, "t": 50, "b": 0, "l": 0}, "scene": {"yaxis": {"title": "", "showgrid": false, "zeroline": false, "showbackground": false, "showticklabels": false, "showline": false}, "xaxis": {"title": "", "showgrid": false, "zeroline": false, "showbackground": false, "showticklabels": false, "showline": false}, "zaxis": {"title": "", "showgrid": false, "zeroline": false, "showbackground": false, "showticklabels": false, "showline": false}}, "hovermode": "closest", "height": 600}, {"linkText": "Export to plot.ly", "showLink": true})});</script>


# Plot specific community as defined by name of boundaries


```python
comm_name = "chr1:12854095-13129135"
pick_comm_graph = dna_int_graph.subgraph(boundary2nodes[comm_name])
pick_comm_graph.graph["bounds"] = comm_name
plot_comm(pick_comm_graph)
```


<div id="4f922876-5656-49a9-93cd-3a0dfa8c8f66" style="height: 600px; width: 600px;" class="plotly-graph-div"></div><script type="text/javascript">require(["plotly"], function(Plotly) { window.PLOTLYENV=window.PLOTLYENV || {};window.PLOTLYENV.BASE_URL="https://plot.ly";Plotly.newPlot("4f922876-5656-49a9-93cd-3a0dfa8c8f66", [{"mode": "lines", "hoverinfo": "none", "line": {"width": 1}, "x": [0.6469172013752597, 0.4326766976156677, null, 0.6469172013752597, 0.8589889131166288, null, 0.3669123367277949, 0.4326766976156677, null, 0.31437632907029095, 0.19973793766864764, null, 0.31437632907029095, 0.30442991513437223, null, 0.31437632907029095, 0.4323488612798084, null, 0.03186776506318509, 0.30442991513437223, null, 0.4326766976156677, 0.30442991513437223, null, 0.4326766976156677, 0.7727139190014364, null, 0.0, 0.30442991513437223, null], "name": "", "type": "scatter", "y": [0.8364318676201977, 0.6378102620377388, null, 0.8364318676201977, 1.0, null, 0.8880351632709166, 0.6378102620377388, null, 0.3190562766283853, 0.0, null, 0.3190562766283853, 0.5178487681043418, null, 0.3190562766283853, 0.1149368940261165, null, 0.692934271995325, 0.5178487681043418, null, 0.6378102620377388, 0.5178487681043418, null, 0.6378102620377388, 0.539683973165994, null, 0.43305716566400326, 0.5178487681043418, null]}, {"mode": "markers", "hoverinfo": "text", "text": ["chr1:13121498-13123090", "chr1:13115518-13119048", "chr1:13093130-13095652", "chr1:13039584-13041392", "chr1:13049045-13050571", "chr1:13127929-13129134", "chr1:13041583-13043936", "chr1:12979738-12983484", "chr1:12854095-12856319", "chr1:13123202-13125339", "chr1:13107646-13108403"], "x": [0.6469172013752597, 0.8589889131166288, 0.3669123367277949, 0.31437632907029095, 0.03186776506318509, 0.19973793766864764, 0.4326766976156677, 0.0, 0.4323488612798084, 0.30442991513437223, 0.7727139190014364], "name": "", "type": "scatter", "marker": {"size": [20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20], "opacity": [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1], "color": ["grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey"]}, "y": [0.8364318676201977, 1.0, 0.8880351632709166, 0.3190562766283853, 0.692934271995325, 0.0, 0.6378102620377388, 0.43305716566400326, 0.1149368940261165, 0.5178487681043418, 0.539683973165994]}], {"title": "chr1:12854095-13129135", "width": 600, "yaxis": {"title": "", "showgrid": false, "zeroline": false, "showbackground": false, "showticklabels": false, "showline": false}, "xaxis": {"title": "", "showgrid": false, "zeroline": false, "showbackground": false, "showticklabels": false, "showline": false}, "showlegend": false, "margin": {"r": 0, "t": 50, "b": 0, "l": 0}, "scene": {"yaxis": {"title": "", "showgrid": false, "zeroline": false, "showbackground": false, "showticklabels": false, "showline": false}, "xaxis": {"title": "", "showgrid": false, "zeroline": false, "showbackground": false, "showticklabels": false, "showline": false}, "zaxis": {"title": "", "showgrid": false, "zeroline": false, "showbackground": false, "showticklabels": false, "showline": false}}, "hovermode": "closest", "height": 600}, {"linkText": "Export to plot.ly", "showLink": true})});</script>


# Print top 10 communities with most nodes and plot #1


```python
comm_sort_ix = np.argsort(bound_list_num_nodes)[::-1]
top_10_comm = [bound_list_names[comm] for comm in comm_sort_ix[0:9]]
for comm in top_10_comm:
    cur_nodes = boundary2nodes[comm]
    bound_split = re.split(r"[-:]",comm)
    print("Community name: " 
              + comm + "\tlinear_length: " + 
              str((int(bound_split[2])-int(bound_split[1]))/1000) + "Kb" +
              "\tnum_nodes: " + str(len(cur_nodes)))

comm_top_graph = dna_int_graph.subgraph(boundary2nodes[top_10_comm[0]])
comm_top_graph.graph["bounds"] = top_10_comm[0]
plot_comm(comm_top_graph)
```

    Community name: chr12:88055955-88394658	linear_length: 338.703Kb	num_nodes: 30
    Community name: chr5:64679876-65270221	linear_length: 590.345Kb	num_nodes: 26
    Community name: chr12:12329799-13222601	linear_length: 892.802Kb	num_nodes: 24
    Community name: chr8:122893930-123142836	linear_length: 248.906Kb	num_nodes: 23
    Community name: chr11:96053696-96484946	linear_length: 431.25Kb	num_nodes: 22
    Community name: chr16:29628679-30237265	linear_length: 608.586Kb	num_nodes: 20
    Community name: chr7:148013968-148362843	linear_length: 348.875Kb	num_nodes: 19
    Community name: chr3:83750397-84347196	linear_length: 596.799Kb	num_nodes: 18
    Community name: chr17:29405078-29652228	linear_length: 247.15Kb	num_nodes: 18



<div id="f9eee291-3d44-4a76-8368-c8fc76d3ab77" style="height: 600px; width: 600px;" class="plotly-graph-div"></div><script type="text/javascript">require(["plotly"], function(Plotly) { window.PLOTLYENV=window.PLOTLYENV || {};window.PLOTLYENV.BASE_URL="https://plot.ly";Plotly.newPlot("f9eee291-3d44-4a76-8368-c8fc76d3ab77", [{"mode": "lines", "hoverinfo": "none", "line": {"width": 1}, "x": [0.4029237905134681, 0.6023809292954053, null, 0.4029237905134681, 0.2235122127677524, null, 0.4029237905134681, 0.5731812483875964, null, 0.4029237905134681, 0.2847244363088823, null, 0.4029237905134681, 0.6653838105962229, null, 0.4002952415160389, 0.5731812483875964, null, 0.4002952415160389, 0.4871303926537024, null, 0.4002952415160389, 0.6232715899270614, null, 0.7002564082641586, 0.6587268923336875, null, 0.2847244363088823, 0.5087921797435132, null, 0.2847244363088823, 0.07723748667862952, null, 0.2847244363088823, 0.0, null, 0.2847244363088823, 0.6587268923336875, null, 0.5087921797435132, 0.6023809292954053, null, 0.5087921797435132, 0.6587268923336875, null, 0.5087921797435132, 0.5731812483875964, null, 0.5087921797435132, 0.6653838105962229, null, 0.5087921797435132, 0.6232715899270614, null, 0.4871303926537024, 0.6587268923336875, null, 0.4871303926537024, 0.6232715899270614, null, 0.44721733011381476, 0.5731812483875964, null, 0.44721733011381476, 0.5412942186108427, null, 0.44721733011381476, 0.6587268923336875, null, 0.7596728391162797, 0.6587268923336875, null, 0.7596728391162797, 0.828917820546032, null, 0.7596728391162797, 0.6232715899270614, null, 0.8921960934920347, 1.0, null, 0.8921960934920347, 0.6653838105962229, null, 0.8921960934920347, 0.6232715899270614, null, 1.0, 0.6587268923336875, null, 0.828917820546032, 0.6587268923336875, null, 0.828917820546032, 0.8600620882252519, null, 0.828917820546032, 0.7234902247699853, null, 0.828917820546032, 0.796312146083792, null, 0.6023809292954053, 0.6587268923336875, null, 0.6023809292954053, 0.6653838105962229, null, 0.6023809292954053, 0.5901051927399336, null, 0.6023809292954053, 0.7821947119575972, null, 0.6023809292954053, 0.5412942186108427, null, 0.6023809292954053, 0.5731812483875964, null, 0.5048261510132174, 0.5731812483875964, null, 0.5048261510132174, 0.6232715899270614, null, 0.8600620882252519, 0.6587268923336875, null, 0.8600620882252519, 0.7168387023107785, null, 0.6653838105962229, 0.6587268923336875, null, 0.6653838105962229, 0.7821947119575972, null, 0.6653838105962229, 0.8143778648360722, null, 0.6653838105962229, 0.5731812483875964, null, 0.6653838105962229, 0.6232715899270614, null, 0.5901051927399336, 0.6587268923336875, null, 0.5901051927399336, 0.6232715899270614, null, 0.8143778648360722, 0.6587268923336875, null, 0.5731812483875964, 0.7234902247699853, null, 0.5731812483875964, 0.5412942186108427, null, 0.5731812483875964, 0.6232715899270614, null, 0.5731812483875964, 0.7168387023107785, null, 0.5731812483875964, 0.6587268923336875, null, 0.5412942186108427, 0.6587268923336875, null, 0.5412942186108427, 0.6232715899270614, null, 0.7145279708052729, 0.6587268923336875, null, 0.7234902247699853, 0.6587268923336875, null, 0.6232715899270614, 0.7168387023107785, null, 0.6232715899270614, 0.6587268923336875, null, 0.6041214480347747, 0.6587268923336875, null, 0.7168387023107785, 0.6587268923336875, null, 0.6587268923336875, 0.796312146083792, null], "name": "", "type": "scatter", "y": [0.5016239404953665, 0.4243331615399618, null, 0.5016239404953665, 0.6849399171244392, null, 0.5016239404953665, 0.3298679635376983, null, 0.5016239404953665, 0.40869407195343604, null, 0.5016239404953665, 0.4531585870121174, null, 0.3289321553662846, 0.3298679635376983, null, 0.3289321553662846, 0.49327858743700387, null, 0.3289321553662846, 0.2975221844172718, null, 0.6484071748323244, 0.3599423827731215, null, 0.40869407195343604, 0.39007331966621595, null, 0.40869407195343604, 0.34327653277405396, null, 0.40869407195343604, 0.46728956395763716, null, 0.40869407195343604, 0.3599423827731215, null, 0.39007331966621595, 0.4243331615399618, null, 0.39007331966621595, 0.3599423827731215, null, 0.39007331966621595, 0.3298679635376983, null, 0.39007331966621595, 0.4531585870121174, null, 0.39007331966621595, 0.2975221844172718, null, 0.49327858743700387, 0.3599423827731215, null, 0.49327858743700387, 0.2975221844172718, null, 0.23655007977565964, 0.3298679635376983, null, 0.23655007977565964, 0.3579029034648088, null, 0.23655007977565964, 0.3599423827731215, null, 0.14264084171789093, 0.3599423827731215, null, 0.14264084171789093, 0.2605134653767809, null, 0.14264084171789093, 0.2975221844172718, null, 0.4286644710603474, 0.4105501630576581, null, 0.4286644710603474, 0.4531585870121174, null, 0.4286644710603474, 0.2975221844172718, null, 0.4105501630576581, 0.3599423827731215, null, 0.2605134653767809, 0.3599423827731215, null, 0.2605134653767809, 0.1620616082646487, null, 0.2605134653767809, 0.27404439248680434, null, 0.2605134653767809, 0.3713088328414199, null, 0.4243331615399618, 0.3599423827731215, null, 0.4243331615399618, 0.4531585870121174, null, 0.4243331615399618, 0.18483505304481082, null, 0.4243331615399618, 0.6444198635303455, null, 0.4243331615399618, 0.3579029034648088, null, 0.4243331615399618, 0.3298679635376983, null, 0.14154572767288554, 0.3298679635376983, null, 0.14154572767288554, 0.2975221844172718, null, 0.1620616082646487, 0.3599423827731215, null, 0.1620616082646487, 0.2044877467000591, null, 0.4531585870121174, 0.3599423827731215, null, 0.4531585870121174, 0.6444198635303455, null, 0.4531585870121174, 0.530088694802342, null, 0.4531585870121174, 0.3298679635376983, null, 0.4531585870121174, 0.2975221844172718, null, 0.18483505304481082, 0.3599423827731215, null, 0.18483505304481082, 0.2975221844172718, null, 0.530088694802342, 0.3599423827731215, null, 0.3298679635376983, 0.27404439248680434, null, 0.3298679635376983, 0.3579029034648088, null, 0.3298679635376983, 0.2975221844172718, null, 0.3298679635376983, 0.2044877467000591, null, 0.3298679635376983, 0.3599423827731215, null, 0.3579029034648088, 0.3599423827731215, null, 0.3579029034648088, 0.2975221844172718, null, 0.0, 0.3599423827731215, null, 0.27404439248680434, 0.3599423827731215, null, 0.2975221844172718, 0.2044877467000591, null, 0.2975221844172718, 0.3599423827731215, null, 0.7055546375171862, 0.3599423827731215, null, 0.2044877467000591, 0.3599423827731215, null, 0.3599423827731215, 0.3713088328414199, null]}, {"mode": "markers", "hoverinfo": "text", "text": ["chr12:88108405-88110967", "chr12:88263008-88265918", "chr12:88211553-88214163", "chr12:88339266-88342495", "chr12:88111152-88114342", "chr12:88170115-88172016", "chr12:88149368-88152459", "chr12:88317104-88318711", "chr12:88357599-88358908", "chr12:88207827-88210549", "chr12:88271786-88274490", "chr12:88314790-88316531", "chr12:88164456-88167041", "chr12:88176811-88178411", "chr12:88259084-88261396", "chr12:88156518-88159055", "chr12:88393253-88394657", "chr12:88285467-88288884", "chr12:88252371-88256201", "chr12:88281448-88284855", "chr12:88244152-88248387", "chr12:88172127-88175257", "chr12:88055955-88057122", "chr12:88225286-88227519", "chr12:88228457-88241520", "chr12:88116117-88117824", "chr12:88144828-88148685", "chr12:88134365-88144307", "chr12:88119959-88121410", "chr12:88218915-88221862"], "x": [0.4029237905134681, 0.6023809292954053, 0.5412942186108427, 0.2847244363088823, 0.5087921797435132, 0.4871303926537024, 0.796312146083792, 0.2235122127677524, 0.0, 0.7596728391162797, 0.8921960934920347, 0.8143778648360722, 0.828917820546032, 0.4002952415160389, 0.5048261510132174, 0.8600620882252519, 0.07723748667862952, 0.6653838105962229, 0.5901051927399336, 1.0, 0.5731812483875964, 0.7002564082641586, 0.7145279708052729, 0.7234902247699853, 0.6232715899270614, 0.6041214480347747, 0.7168387023107785, 0.6587268923336875, 0.7821947119575972, 0.44721733011381476], "name": "", "type": "scatter", "marker": {"size": [20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20], "opacity": [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1], "color": ["grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey"]}, "y": [0.5016239404953665, 0.4243331615399618, 0.3579029034648088, 0.40869407195343604, 0.39007331966621595, 0.49327858743700387, 0.3713088328414199, 0.6849399171244392, 0.46728956395763716, 0.14264084171789093, 0.4286644710603474, 0.530088694802342, 0.2605134653767809, 0.3289321553662846, 0.14154572767288554, 0.1620616082646487, 0.34327653277405396, 0.4531585870121174, 0.18483505304481082, 0.4105501630576581, 0.3298679635376983, 0.6484071748323244, 0.0, 0.27404439248680434, 0.2975221844172718, 0.7055546375171862, 0.2044877467000591, 0.3599423827731215, 0.6444198635303455, 0.23655007977565964]}], {"title": "chr12:88055955-88394658", "width": 600, "yaxis": {"title": "", "showgrid": false, "zeroline": false, "showbackground": false, "showticklabels": false, "showline": false}, "xaxis": {"title": "", "showgrid": false, "zeroline": false, "showbackground": false, "showticklabels": false, "showline": false}, "showlegend": false, "margin": {"r": 0, "t": 50, "b": 0, "l": 0}, "scene": {"yaxis": {"title": "", "showgrid": false, "zeroline": false, "showbackground": false, "showticklabels": false, "showline": false}, "xaxis": {"title": "", "showgrid": false, "zeroline": false, "showbackground": false, "showticklabels": false, "showline": false}, "zaxis": {"title": "", "showgrid": false, "zeroline": false, "showbackground": false, "showticklabels": false, "showline": false}}, "hovermode": "closest", "height": 600}, {"linkText": "Export to plot.ly", "showLink": true})});</script>

