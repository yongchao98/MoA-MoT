import pandas as pd
from neuprint import Client

def find_main_asj_target():
    """
    Finds the main projection target of ASJ axons in the C. elegans
    adult male connectome by synapse number.
    """
    try:
        # It is recommended to have your neuprint token stored in a file
        # and the path to that file set as the NEUPRINT_APPLICATION_CREDENTIALS
        # environment variable.
        # See: https://connectome-neuprint.github.io/neuprint-python/docs/quickstart.html
        c = Client('neuprint.janelia.org', dataset='c-elegans-male')
        c.fetch_version()
    except Exception as e:
        print(f"Failed to connect to neuPrint server. Please ensure you have an internet connection and a valid token.")
        print(f"Error: {e}")
        return

    # Define the presynaptic neurons of interest (ASJ left and right)
    presynaptic_neurons = ['ASJL', 'ASJR']

    # Fetch all outgoing connections from ASJ neurons
    try:
        connections_df = c.fetch_simple_connections(pre_criteria={'name': presynaptic_neurons})
        if connections_df is None or connections_df.empty:
            print(f"No connections found for neurons: {presynaptic_neurons}")
            return
    except Exception as e:
        print(f"An error occurred while fetching data: {e}")
        return

    # Group by the postsynaptic neuron name and sum the synapse counts (weight)
    target_summary = connections_df.groupby('post_name')['weight'].sum().reset_index()

    # Sort to find the main target with the most synapses
    target_summary_sorted = target_summary.sort_values(by='weight', ascending=False)

    # Get the top target cell
    if not target_summary_sorted.empty:
        main_target_info = target_summary_sorted.iloc[0]
        target_cell = main_target_info['post_name']
        total_synapses = main_target_info['weight']

        # Find the breakdown of synapses from ASJL and ASJR to the main target
        target_connections = connections_df[connections_df['post_name'] == target_cell]
        
        synapses_from_asjl = 0
        synapses_from_asjr = 0

        # There might be no connection from one of the pair
        if not target_connections[target_connections['pre_name'] == 'ASJL'].empty:
            synapses_from_asjl = target_connections[target_connections['pre_name'] == 'ASJL']['weight'].iloc[0]
        
        if not target_connections[target_connections['pre_name'] == 'ASJR'].empty:
            synapses_from_asjr = target_connections[target_connections['pre_name'] == 'ASJR']['weight'].iloc[0]

        print(f"The main projection target of ASJ axons is the cell: {target_cell}")
        print("\nThis is based on the total number of synapses received from both ASJL and ASJR.")
        print("\nEquation:")
        print(f"{synapses_from_asjl} (from ASJL) + {synapses_from_asjr} (from ASJR) = {total_synapses} synapses to {target_cell}")
    else:
        print("Could not determine the main target.")

if __name__ == '__main__':
    find_main_asj_target()