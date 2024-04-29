import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# Load data and matrices
def CreateMatrices(filename='CRT.xlsx'):
    df = pd.read_excel(filename)
    NodeNames = df['Nodes'].tolist()
    NumOfNodes = len(NodeNames)
    Mact = np.zeros((NumOfNodes, NumOfNodes))
    Minh = np.zeros((NumOfNodes, NumOfNodes))
    Clamped = np.zeros(NumOfNodes)

    for i, row in df.iterrows():
        activators = str(row['Activators']).split(',')
        inhibitors = str(row['Inhibitors']).split(',')

        for act in activators:
            if act in NodeNames:
                Mact[i, NodeNames.index(act)] = 1

        for inh in inhibitors:
            if inh in NodeNames:
                Minh[i, NodeNames.index(inh)] = 1

    return Mact, Minh, NodeNames, Clamped

# Load data and matrices
Mact, Minh, NodeNames, Clamped = CreateMatrices()

# Create a directed graph
G = nx.DiGraph()

# Add nodes to the graph
for name in NodeNames:
    G.add_node(name)

# Add edges based on activation and inhibition matrices
for i, row in enumerate(Mact):
    for j, val in enumerate(row):
        if val == 1:
            G.add_edge(NodeNames[i], NodeNames[j], color='blue')  # Activator edge
for i, row in enumerate(Minh):
    for j, val in enumerate(row):
        if val == 1:
            G.add_edge(NodeNames[j], NodeNames[i], color='red')  # Inhibitor edge

# Draw the network
pos = nx.spring_layout(G, k=0.5, iterations=50, seed=42)  # Adjust layout parameters
plt.figure(figsize=(12, 8))  # Adjust figure size
nx.draw(G, pos, with_labels=True, node_color='lightblue', node_size=1500, font_size=10, font_weight='bold')
edge_colors = nx.get_edge_attributes(G, 'color').values()
nx.draw_networkx_edges(G, pos, edge_color=edge_colors)

# Save the network image
plt.savefig("network_image.png", dpi=300)  # Increase dpi for higher resolution
plt.show()
