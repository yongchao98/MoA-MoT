import pandas as pd
try:
    from celegans_metar.connectome import get_connectome
except ImportError:
    print("Please install the required library by running: pip install celegans-metar-connectome")
    exit()

def find_main_synaptic_target():
    """
    Finds the main synaptic target of a specified neuron in the C. elegans connectome.
    """
    try:
        # 1. Load the adult hermaphrodite connectome data from Witvliet et al., 2021
        connectome = get_connectome()
    except Exception as e:
        print(f"Failed to download or load connectome data. Please check your internet connection.")
        print(f"Error: {e}")
        return

    # 2. Define the presynaptic neurons. The question asks about ASJ, which is a pair (ASJL/ASJR).
    pre_neurons = ['ASJL', 'ASJR']

    # 3. Filter for connections originating from ASJ neurons.
    asj_outputs = connectome[connectome['pre_cell'].isin(pre_neurons)]

    # 4. Group by the postsynaptic partner and sum the number of synapses ('N').
    target_synapse_counts = asj_outputs.groupby('post_cell')['N'].sum()

    # 5. Sort the results to find the main target.
    sorted_targets = target_synapse_counts.sort_values(ascending=False)

    if sorted_targets.empty:
        print(f"No synaptic targets found for {pre_neurons}.")
        return

    # 6. Get the top target neuron and the number of synapses.
    main_target_cell = sorted_targets.index[0]
    main_target_synapses = int(sorted_targets.iloc[0]) # Convert to int for clean printing

    # 7. Print the results, showing the contribution from ASJ to the top target.
    print(f"Finding the main projection target of ASJ (ASJL + ASJR) axons by synapse number.")
    print("-" * 50)
    print(f"Total synapses from ASJ (ASJL + ASJR) to {main_target_cell}: {main_target_synapses}")
    print("-" * 50)
    print(f"The main projection target of ASJ axons is the neuron '{main_target_cell}'.")


if __name__ == "__main__":
    find_main_synaptic_target()
