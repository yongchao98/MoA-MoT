import pandas as pd
try:
    import navis
except ImportError:
    print("Please install the 'navis' library by running: pip install navis")
    exit()

def find_main_asj_target():
    """
    Finds the main synaptic target of ASJ neurons in the adult C. elegans connectome.
    """
    try:
        # Fetch the Witvliet et al. 2021 adult connectome
        # This dataset contains synaptic connectivity data.
        conn = navis.fetch_connectome('witvliet-2021')
    except Exception as e:
        print(f"Failed to fetch connectome data. Please check your internet connection.")
        print(f"Error: {e}")
        return

    # Define the source neurons (the ASJ pair)
    source_neurons = ['ASJL', 'ASJR']

    # Filter the connectome to get all outputs from ASJ neurons
    # 'pre_name' is the presynaptic (source) neuron
    asj_outputs = conn[conn['pre_name'].isin(source_neurons)]

    if asj_outputs.empty:
        print("Could not find any outputs from ASJ neurons in the dataset.")
        return

    # Group by the postsynaptic (target) neuron and sum the synapse counts ('N')
    # This gives us the total number of synapses from ASJ to each target cell.
    target_summary = asj_outputs.groupby('post_name')['N'].sum()

    # Find the target cell with the maximum number of synapses
    main_target_cell = target_summary.idxmax()
    total_synapses = target_summary.max()

    # Get the individual synapse counts from ASJL and ASJR to the main target
    # We query the filtered dataframe for the specific connections.
    # .sum() is used in case there are multiple entries for the same pair, though unlikely.
    asjl_to_target_synapses = asj_outputs[
        (asj_outputs['pre_name'] == 'ASJL') & (asj_outputs['post_name'] == main_target_cell)
    ]['N'].sum()

    asjr_to_target_synapses = asj_outputs[
        (asj_outputs['pre_name'] == 'ASJR') & (asj_outputs['post_name'] == main_target_cell)
    ]['N'].sum()

    # Print the final result and the breakdown
    print(f"The main projection target of ASJ axons is the {main_target_cell} neuron.")
    print("\nThe synaptic equation is:")
    print(f"ASJL -> {main_target_cell} = {asjl_to_target_synapses}")
    print(f"ASJR -> {main_target_cell} = {asjr_to_target_synapses}")
    print("--------------------")
    print(f"ASJ (Total) -> {main_target_cell} = {total_synapses}")

if __name__ == "__main__":
    find_main_asj_target()