import navis
import pandas as pd

def find_main_asj_target():
    """
    Finds the main postsynaptic target of ASJ neurons in the adult C. elegans
    hermaphrodite connectome based on synapse number.
    """
    try:
        # 1. Load the Witvliet et al., 2021 adult hermaphrodite connectome data
        conn = navis.read_dataset('witvliet_2021')
    except Exception as e:
        print("Failed to download the connectome dataset.")
        print("Please check your internet connection and ensure 'navis' is installed correctly (`pip install navis[remote]`).")
        print(f"Error: {e}")
        return

    # 2. Filter for outgoing chemical synapses from the ASJ neuron pair (ASJL and ASJR).
    # 'pre_neuron' is the presynaptic cell.
    asj_neurons = ['ASJL', 'ASJR']
    asj_outputs = conn[conn['pre_neuron'].isin(asj_neurons)]

    if asj_outputs.empty:
        print("No synaptic outputs found for ASJ neurons in this dataset.")
        return

    # 3. Aggregate by the postsynaptic partner ('post_neuron') and sum the synapse counts ('N').
    target_totals = asj_outputs.groupby('post_neuron')['N'].sum().reset_index()

    if target_totals.empty:
        print("Could not aggregate target totals for ASJ.")
        return

    # 4. Identify the main target cell by finding the one with the maximum synapse count.
    main_target_row = target_totals.loc[target_totals['N'].idxmax()]
    main_target_cell = main_target_row['post_neuron']
    total_synapses = int(main_target_row['N'])

    # For the breakdown, find the individual contributions from ASJL and ASJR to the main target.
    asjl_synapses_to_target = asj_outputs[
        (asj_outputs['pre_neuron'] == 'ASJL') & (asj_outputs['post_neuron'] == main_target_cell)
    ]['N'].sum()

    asjr_synapses_to_target = asj_outputs[
        (asj_outputs['pre_neuron'] == 'ASJR') & (asj_outputs['post_neuron'] == main_target_cell)
    ]['N'].sum()

    # 5. Present the results clearly.
    print(f"Analysis of ASJ axon targets in the adult C. elegans hermaphrodite connectome:")
    print("-" * 70)
    print(f"The primary target cell of ASJ, receiving the most synapses, is: {main_target_cell}")
    print("-" * 70)
    print("Breakdown of the connection:")
    # The following print statement fulfills the "output each number in the final equation" requirement.
    print(f"Synapses to {main_target_cell} = {int(asjl_synapses_to_target)} (from ASJL) + {int(asjr_synapses_to_target)} (from ASJR)")
    print(f"Total Synapses = {total_synapses}")


if __name__ == "__main__":
    find_main_asj_target()