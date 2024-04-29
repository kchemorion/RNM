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
        # Continue with simulation setup and plotting...
    except Exception as e:
        print("Failed to load the model:", str(e))
    
    r.reset()  # Reset the model to initial conditions

    # Define simulation parameters
    start_time = 0
    end_time = 30
    num_points = 100

    # Simulate the model
    results = r.simulate(start_time, end_time, num_points)

    # Correct indexing - MATLAB to Python (MATLAB 1-based, Python 0-based)
    # Define node groups for plotting based on earlier MATLAB indices adjusted for Python
    pro_inflammatory = ['IL6', 'TNFa', 'IL1b']  # Example names, check your actual node names
    anti_inflammatory = ['IL4', 'IL10']
    growth_factors = ['BMP2', 'FGF2', 'TGFB1']
    ecm_destruction = ['ADAMTS4', 'MMP1', 'MMP3', 'MMP13']
    ecm_synthesis = ['ACAN', 'COL2A1']
    hypertrophy = ['COL10A1', 'MMP13']  # Placeholder names

    # Plot results
    plt.figure(figsize=(12, 8))

    def plot_group(index, group, title, filename):
        plt.subplot(2, 3, index)
        for node in group:
            if node in results.colnames:
                plt.plot(results['time'], results[node], label=node)
        plt.title(title)
        plt.xlabel('Time')
        plt.ylabel('Concentration')
        plt.legend()
        plt.savefig(filename)  # Save each plot to a file

    plot_group(1, pro_inflammatory, 'Pro-Inflammatory', 'pro_inflammatory.png')
    plot_group(2, anti_inflammatory, 'Anti-Inflammatory', 'anti_inflammatory.png')
    plot_group(3, growth_factors, 'Growth Factors', 'growth_factors.png')
    plot_group(4, ecm_destruction, 'ECM Destruction', 'ecm_destruction.png')
    plot_group(5, ecm_synthesis, 'ECM Synthesis', 'ecm_synthesis.png')
    plot_group(6, hypertrophy, 'Hypertrophy', 'hypertrophy.png')

    plt.tight_layout()
    plt.show()

    # Export the model to SBML
    sbml_model = r.getSBML()
    sbml_path = 'model.sbml'
    with open(sbml_path, 'w') as sbml_file:
        sbml_file.write(sbml_model)

    print(f"SBML model exported successfully to {sbml_path}.")
