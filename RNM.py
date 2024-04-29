import tellurium as te
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Function to Create Matrices from Excel
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

# Function to create an Antimony model
def sanitize_name(name):
    # Replace special characters with underscore
    return name.replace("-", "_").replace("/", "_").replace("β", "beta").replace("γ", "gamma").replace("α", "alpha")

def createAntimonyModel(Mact, Minh, NodeNames, Clamped):
    model_str = 'model networkODE()\n'
    
    # Sanitize and define nodes as species in the model
    NodeNames = [sanitize_name(name) for name in NodeNames]
    
    for name in NodeNames:
        model_str += f'    var {name};\n'  # Declare variables
        model_str += f'    {name} = 1;  // Initial condition\n'
    
    model_str += '\n'
    
    # Add ODEs for each node based on activators and inhibitors
    for i, node in enumerate(NodeNames):
        activator_indices = np.where(Mact[i] == 1)[0]
        inhibitor_indices = np.where(Minh[i] == 1)[0]

        activation_str = ' + '.join([f'{NodeNames[j]}' for j in activator_indices])
        inhibition_str = ' + '.join([f'(1 - {NodeNames[j]} / (1 + {NodeNames[j]}))' for j in inhibitor_indices])

        if activator_indices.size > 0:
            activation_str = f'({activation_str}) / (1 + {activation_str})'
        else:
            activation_str = '0'

        if inhibitor_indices.size > 0:
            inhibition_str = f'({inhibition_str})'
        else:
            inhibition_str = '1'

        model_str += f'    {node}\' = {activation_str} * {inhibition_str} - {node};\n'

    model_str += 'end'
    
    return model_str

# Main execution block
if __name__ == "__main__":
    # Load data and matrices
    Mact, Minh, NodeNames, Clamped = CreateMatrices()

    # Generate the Antimony model string
    antimony_model = createAntimonyModel(Mact, Minh, NodeNames, Clamped)

    print("Generated Antimony Model:")
    print(antimony_model)

    # Load the model in roadrunner
    try:
        r = te.loada(antimony_model)
        print("Model loaded successfully.")
    except Exception as e:
        print("Failed to load the model:", str(e))
    
    r.reset()  # Reset the model to initial conditions

    # Define simulation parameters
    start_time = 0
    end_time = 30
    num_points = 100

    # Simulate the model
    results = r.simulate(start_time, end_time, num_points)

    # Node group definitions - assuming indices are accurate and NodeNames are correctly ordered
    pro_inflammatory = [NodeNames[i] for i in [0, 1, 2]]  # Adjusted indices to Python's 0-based indexing
    anti_inflammatory = [NodeNames[i] for i in [3]]
    growth_factors = [NodeNames[i] for i in [9]]
    ecm_destruction = [NodeNames[i] for i in [0, 1]]
    ecm_synthesis = [NodeNames[i] for i in [2]]
    hypertrophy = [NodeNames[i] for i in [49]]  # Adjust index for Python

    # Plotting results
    plt.figure(figsize=(18, 10))

    # Helper function to plot each group
    def plot_group(data, group, subplot_index, title):
        plt.subplot(2, 3, subplot_index)
        for node in group:
            if node in data.colnames:
                plt.plot(data['time'], data[node], label=node)
        plt.title(title)
        plt.xlabel('Time')
        plt.ylabel('Normalized Presence')
        plt.legend()
        plt.savefig(f"{title.replace(' ', '_').lower()}.png")  # Save figure to file

    plot_group(results, pro_inflammatory, 1, 'Pro-Inflammatory')
    plot_group(results, anti_inflammatory, 2, 'Anti-Inflammatory')
    plot_group(results, growth_factors, 3, 'Growth Factors')
    plot_group(results, ecm_destruction, 4, 'ECM Destruction')
    plot_group(results, ecm_synthesis, 5, 'ECM Synthesis')
    plot_group(results, hypertrophy, 6, 'Hypertrophy')

    plt.tight_layout()
    plt.show()

    # Export the model to SBML
    sbml_model = r.getSBML()
    with open('model.sbml', 'w') as sbml_file:
        sbml_file.write(sbml_model)
    print("SBML model exported and plots saved successfully.")

