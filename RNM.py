import numpy as np
import tellurium as te
import matplotlib.pyplot as plt
import os
import libsbml
import pygraphviz as pgv
import pandas as pd

def ODESysFun(t, X):
    # Parameter Definition
    NumOfNodes = len(X)
    gamma = np.ones(NumOfNodes)
    h = 10  # Steepness of activation
    Mact, Minh, Clamped = CreateMatrices()

    w = np.zeros(NumOfNodes)
    f = np.zeros(NumOfNodes)

    # Calculation of ODEs from activation and inhibition matrices
    for i in range(NumOfNodes):
        Ract = Mact[i]
        Rinh = Minh[i]

        if (np.any(Rinh) == 0 and np.any(Ract) == 1):  # Node i has no inhibitors and at least 1 activator
            sum_alpha = np.sum(Ract)
            sum_alpha_X = np.dot(Ract, X)

            w[i] = ((1 + sum_alpha) / sum_alpha) * (sum_alpha_X / (1 + sum_alpha_X))

        elif (np.any(Ract) == 0 and np.any(Rinh) == 1):  # Node i has no activators and at least 1 inhibitor
            sum_beta = np.sum(Rinh)
            sum_beta_X = np.dot(Rinh, X)

            w[i] = 1 - ((1 + sum_beta) / sum_beta) * (sum_beta_X / (1 + sum_beta_X))

        elif (np.any(Ract) == 1 and np.any(Rinh) == 1):  # Node i has inhibitors and activators
            sum_alpha = np.sum(Ract)
            sum_beta = np.sum(Rinh)

            sum_alpha_X = np.dot(Ract, X)
            sum_beta_X = np.dot(Rinh, X)

            w[i] = ((1 + sum_alpha) / sum_alpha) * (sum_alpha_X / (1 + sum_alpha_X)) * (
                        1 - ((1 + sum_beta) / sum_beta) * (sum_beta_X / (1 + sum_beta_X)))
        else:
            w[i] = 0

        f[i] = (-np.exp(0.5 * h) + np.exp(-h * (w[i] - 0.5))) / ((1 - np.exp(0.5 * h)) * (1 + np.exp(-h * (w[i] - 0.5)))) - gamma[i] * X[i]
        if Clamped[i] == 1:
            f[i] = 0

    return f

def CreateMatrices():
    # Read data from CRT.xlsx file
    df = pd.read_excel('CRT.xlsx')

    NodeNames = df['Nodes'].tolist()
    NumOfNodes = len(NodeNames)
    Mact = np.zeros((NumOfNodes, NumOfNodes))
    Minh = np.zeros((NumOfNodes, NumOfNodes))
    Clamped = np.zeros(NumOfNodes)

    for i, row in df.iterrows():
        activators = str(row['Activators']).split(',')
        inhibitors = str(row['Inhibitors']).split(',')

        for act in activators:
            if act:
                Mact[i, NodeNames.index(act)] = 1

        for inh in inhibitors:
            if inh:
                Minh[i, NodeNames.index(inh)] = 1

    return Mact, Minh, Clamped

# Solve the ODE system
t = np.linspace(0, 30, 100)  # Time points
Xinit = np.random.rand(61)
r = te.loada('')  # Load your model here
result = r.simulate(0, 30, 100)
Xout = result[:, 1:]  # Get simulation results for species only

# Plot the results
pro_infl = [0, 1, 2]  # Indices of nodes for pro-inflammatory factors
anti_infl = [3]  # Indices of nodes for anti-inflammatory factors
growth = [10]  # Indices of nodes for growth factors
ecm_destr = [0, 1]  # Indices of nodes for ECM destruction factors
ecm_genesis = [2]  # Indices of nodes for ECM synthesis factors
hypertrophy = [50]  # Indices of nodes for hypertrophy

plt.figure(figsize=(14, 10))

print("Simulation Results:")
print(result)

print("Species Concentrations:")
print(Xout)


plt.subplot(2, 3, 1)
plt.plot(result[:, 0], Xout[:, pro_infl])
plt.title('Pro-Inflammatory')
plt.xlabel('Time')
plt.ylabel('Concentration')
plt.legend(['Pro-Inflammatory'])

# Plot other subplots similarly

plt.tight_layout()

# Save the plots
plt.savefig('plot.png')

# Export the model to SBML
r.exportToSBML('network_model.sbml')

# Attempt to draw the network
try:
    sbml_document = libsbml.readSBMLFromFile('network_model.sbml')
    sbml_model = sbml_document.getModel()
    graph = pgv.AGraph(strict=False, directed=True)

    # Populate the graph with nodes (species) and edges (reactions)
    for species in sbml_model.getListOfSpecies():
        graph.add_node(species.getId(), label=species.getName())

    for reaction in sbml_model.getListOfReactions():
        for reactant in reaction.getListOfReactants():
            for product in reaction.getListOfProducts():
                graph.add_edge(reactant.getSpecies(), product.getSpecies())

    # Apply layout and save the graph
    graph.layout(prog='dot')  # 'dot' layout for directed graphs
    graph.draw('network_visualization.png')

    print("Network visualization saved as 'network_visualization.png'.")
except Exception as e:
    print("Failed to draw network graph:", str(e))