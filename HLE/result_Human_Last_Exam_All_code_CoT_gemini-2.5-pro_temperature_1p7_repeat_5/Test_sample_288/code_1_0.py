import pandas as pd
from c_elegans_neuprint.client import Client
from neuprint import fetch_synapse_connections

def find_main_synapse_target():
    """
    Connects to the c-elegans-neuprint database to find the main synaptic
    target of the ASJ neuron based on synapse count.
    """
    try:
        # Connect to the neuprint server with the C. elegans dataset
        # Dataset from: Cook et al., 2019; updated 2023-03-05
        c = Client('https://neuprint.janelia.org', dataset='celegans_20230305')
        
        # Define the presynaptic (source) neuron
        source_neuron = 'ASJ'
        
        # Fetch all outgoing synaptic connections from the ASJ neuron
        connections_df = fetch_synapse_connections(source_criteria={'instance': source_neuron})
        
        # Ensure we found connections and the dataframe is not empty
        if connections_df.empty:
            print(f"No outgoing connections were found for neuron {source_neuron}.")
            return

        # Group connections by the postsynaptic (target) neuron instance and sum the 'weight'
        # The 'weight' column represents the number of synaptic contacts.
        target_summary = connections_df.groupby('instance_post')['weight'].sum().reset_index()
        
        # Sort the results to find the main target
        target_summary_sorted = target_summary.sort_values('weight', ascending=False)
        
        # Get the top target (the one with the most synapses)
        main_target = target_summary_sorted.iloc[0]
        
        target_cell = main_target['instance_post']
        synapse_count = main_target['weight']
        
        print(f"Querying the adult C. elegans connectome for the main target of neuron ASJ...")
        print(f"The primary projection target of ASJ, by synapse number, is {target_cell} with {synapse_count} synapses.")

    except Exception as e:
        print(f"An error occurred while querying the neuPrint database: {e}")
        print("Please check your internet connection and the status of the neuprint.janelia.org server.")

if __name__ == "__main__":
    find_main_synapse_target()