import pandas as pd

def find_main_synapse_target():
    """
    This script finds the main postsynaptic target of ASJ neurons in C. elegans
    based on the number of chemical synapses, using the Cook et al. 2019 connectome data.
    """
    try:
        # URL to the Cook et al., 2019 connectome dataset
        url = "https://raw.githubusercontent.com/c-elegans-open-data/c-elegans-open-data/master/Cook_et_al_2019/suppl_files/SI-2-connectome.csv"

        # Load the dataset into a pandas DataFrame
        df = pd.read_csv(url)

        # Define the presynaptic neurons of interest (ASJ left and right)
        presynaptic_neurons = ['ASJL', 'ASJR']

        # Filter for chemical synapses ('S') from ASJL or ASJR
        asj_axon_outputs = df[
            (df['Neuron 1'].isin(presynaptic_neurons)) &
            (df['Type'] == 'S')
        ]

        if asj_axon_outputs.empty:
            print("Could not find any synaptic outputs for ASJ neurons in the dataset.")
            return

        # Group by the postsynaptic neuron ('Neuron 2') and sum the synapse counts ('N')
        target_synapse_counts = asj_axon_outputs.groupby('Neuron 2')['N'].sum()

        if target_synapse_counts.empty:
            print("No postsynaptic targets found for ASJ neurons.")
            return

        # Find the target cell with the maximum number of synapses
        main_target_cell = target_synapse_counts.idxmax()
        max_synapses = target_synapse_counts.max()

        # Print the final result, showing the target cell and the synapse count
        print(f"The main projection target of ASJ axons is the cell {main_target_cell}.")
        print(f"Total number of synapses from ASJ to {main_target_cell}: {max_synapses}")


    except Exception as e:
        print(f"An error occurred: {e}")
        print("Please check your internet connection or the data URL.")

# Execute the function
find_main_synapse_target()