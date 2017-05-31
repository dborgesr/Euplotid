
## Pipeline to find interconnected regions
This pipeline takes in an <a href="https://github.com/danielsday/origami/"> Origami </a> processed DNA-DNA interactions generated using [fq2ChIAInts](fq2ChIAInts.md) and finds interconnected networks of nodes (communities) using default or picked resolution. Plotting functions are also included to visualize output  
Range for resolution is from [0,1] with 1e-4 being something close to Insulated Neighborhoods, higher resolution makes bigger communities. You can also pick the "level", we start at level 0 then merging level 0 communities we get level 1 communities, so on and so forth.


```python
# take a peek to see what DNA-DNA interactions are available
!ls /input_dir/*.bedpe -lt
#!conda install -y -c conda-forge  lua=5.3.3 
```


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

# Plot specific community as defined by name of boundaries


```python
comm_name = "chr1:12854095-13129135"
pick_comm_graph = dna_int_graph.subgraph(boundary2nodes[comm_name])
pick_comm_graph.graph["bounds"] = comm_name
plot_comm(pick_comm_graph)
```

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
