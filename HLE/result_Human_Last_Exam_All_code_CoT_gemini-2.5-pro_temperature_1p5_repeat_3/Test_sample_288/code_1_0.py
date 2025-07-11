import pandas as pd
from celegans_connectome import connectome

def find_main_target():
    """
    This function identifies the main synaptic target of the ASJ neurons
    in the adult C. elegans connectome.
    """
    try:
        # Load the connectome data for the adult hermaphrodite.
        # This dataset contains synapse information, including pre- and post-synaptic cells
        # and the number of synapses ('N') for each connection.
        synapse_df = connectome.get_data()
    except Exception as e:
        print(f"Could not load the connectome data. Please ensure you have an internet connection")
        print(f"and the 'c-elegans-connectome' library is installed ('pip install c-elegans-connectome').")
        print(f"Error: {e}")
        return

    # Define the source neurons. ASJ is a left-right pair.
    source_neurons = ['ASJL', 'ASJR']

    # Filter the data to get all synaptic outputs from the ASJ neurons.
    asj_outputs = synapse_df[synapse_df['pre_cell'].isin(source_neurons)]

    # Group the filtered connections by the target cell ('post_cell') and
    # sum the number of synapses ('N') from both ASJL and ASJR.
    target_synapse_counts = asj_outputs.groupby('post_cell')['N'].sum()

    # Find the target cell that receives the most synapses from ASJ.
    # idxmax() gives the name of the cell, and max() gives the synapse count.
    main_target_cell = target_synapse_counts.idxmax()
    max_synapses = target_synapse_counts.max()

    # Print the result.
    print(f"The source neuron pair is ASJ (ASJL & ASJR).")
    print(f"The main projection target of ASJ in terms of synapse number is the neuron: {main_target_cell}")
    print(f"ASJ makes a total of {max_synapses} synapses onto {main_target_cell}.")

if __name__ == "__main__":
    find_main_target()
<<<AIA>>>