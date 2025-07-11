import pandas as pd

def find_main_synaptic_target():
    """
    This function downloads the C. elegans connectome data, identifies the main
    synaptic target of ASJ neurons, and prints the result.
    """
    # The URL for the adult hermaphrodite connectome edge list from Cook et al. (2019)
    url = "https://raw.githubusercontent.com/c-elegans-connectome/celegans-connectome-data/master/herm_full_edgelist.csv"

    try:
        # Load the connectome data into a pandas DataFrame
        connectome_df = pd.read_csv(url)
    except Exception as e:
        print(f"Error: Could not download or read the data from the URL.")
        print(f"Details: {e}")
        return

    # Filter for all synaptic outputs from the ASJ neuron pair (ASJL and ASJR)
    asj_outputs = connectome_df[connectome_df['pre'].isin(['ASJL', 'ASJR'])]

    if asj_outputs.empty:
        print("Could not find any synaptic outputs from ASJ neurons in the dataset.")
        return

    # Group by the postsynaptic ('post') neuron and sum the number of synapses ('n')
    target_synapse_counts = asj_outputs.groupby('post')['n'].sum()

    # Find the target neuron that receives the maximum number of synapses
    main_target_neuron = target_synapse_counts.idxmax()
    total_synapses_to_target = target_synapse_counts.max()

    # To show the breakdown, find the individual contributions from ASJL and ASJR
    # to the main target. We use .sum() to handle cases where a connection might not exist.
    synapses_from_asjl = asj_outputs[(asj_outputs['pre'] == 'ASJL') & (asj_outputs['post'] == main_target_neuron)]['n'].sum()
    synapses_from_asjr = asj_outputs[(asj_outputs['pre'] == 'ASJR') & (asj_outputs['post'] == main_target_neuron)]['n'].sum()

    # Print the final results
    print(f"Based on the Cook et al. (2019) dataset, the main projection target of ASJ axons is: {main_target_neuron}")
    print("\n--- Synapse Calculation ---")
    print(f"Synapses from ASJL to {main_target_neuron}: {synapses_from_asjl}")
    print(f"Synapses from ASJR to {main_target_neuron}: {synapses_from_asjr}")
    print(f"Total synapses from ASJ to {main_target_neuron}: {synapses_from_asjl} + {synapses_from_asjr} = {int(total_synapses_to_target)}")

if __name__ == '__main__':
    find_main_synaptic_target()