import cconnectome as cc
import pandas as pd

# Set pandas to display all rows, to ensure the full list of targets is shown if needed
# pandas.set_option('display.max_rows', None)

def find_main_synapse_target():
    """
    Finds the main synaptic target of a given neuron in the C. elegans connectome.
    """
    try:
        # Step 1: Access the connectome data
        print("Connecting to the C. elegans connectome dataset...")
        synapses_df = cc.get_synapses()
        print("Successfully loaded synapse data.")

        # Step 2: Define the source neurons and filter connections
        source_neurons = ['ASJL', 'ASJR']
        print(f"Filtering for all synaptic outputs from {', '.join(source_neurons)}...")
        asj_outputs_df = synapses_df[synapses_df['pre_name'].isin(source_neurons)]

        if asj_outputs_df.empty:
            print(f"No synaptic outputs found for {', '.join(source_neurons)}.")
            return

        # Step 3 & 4: Group by target neuron and sum synapse counts
        print("Aggregating synapse counts for each target neuron...")
        # 'weight' column represents the number of synapses in a connection
        target_synapse_counts = asj_outputs_df.groupby('post_name')['weight'].sum().sort_values(ascending=False)

        # Step 5: Identify the main target
        if target_synapse_counts.empty:
            print(f"Could not determine synapse counts for targets of {', '.join(source_neurons)}.")
            return

        main_target_neuron = target_synapse_counts.index[0]
        max_synapses = target_synapse_counts.iloc[0]

        # Step 6: Print the result
        print("\n--- Result ---")
        print(f"The main projection target of ASJ ({', '.join(source_neurons)}) axons is the neuron '{main_target_neuron}'.")
        print(f"The total number of synapses from ASJ to {main_target_neuron} is {max_synapses}.")

    except Exception as e:
        print(f"An error occurred: {e}")
        print("Please ensure you have the 'cconnectome' library installed ('pip install cconnectome') and a working internet connection.")

if __name__ == '__main__':
    find_main_synapse_target()