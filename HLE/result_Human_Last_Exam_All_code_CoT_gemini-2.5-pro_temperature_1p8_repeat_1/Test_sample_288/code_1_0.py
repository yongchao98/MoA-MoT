import pandas as pd

def find_asj_target():
    """
    Finds the main postsynaptic target of ASJ neurons in C. elegans.
    This function uses a publicly available connectome dataset to calculate
    which neuron receives the most synapses from the ASJ neuron pair (ASJL/ASJR).
    """
    try:
        from celegans_connectome.load import load_connectome
    except ImportError:
        print("Error: The 'celegans-connectome' library is not installed.")
        print("Please install it to run this code:")
        print("pip install celegans-connectome pandas")
        return

    # Load the chemical synapse data for the adult hermaphrodite from Kato et al., 2014
    # The 'chem' DataFrame contains columns: 'pre' (presynaptic), 'post' (postsynaptic), 'synapses'
    try:
        chem, _ = load_connectome(dataset_name="c-elegans-splinter-connectome")
    except Exception as e:
        print(f"Failed to load connectome data. Error: {e}")
        print("Please ensure you have an internet connection to download the data.")
        return

    # Filter for all connections originating from ASJL or ASJR
    asj_axons = ['ASJL', 'ASJR']
    asj_connections = chem[chem['pre'].isin(asj_axons)]

    if asj_connections.empty:
        print("No connections from ASJ neurons were found in the dataset.")
        return

    # Group connections by the target (postsynaptic) neuron and sum the synapses
    target_synapse_counts = asj_connections.groupby('post')['synapses'].sum()

    # Find the target neuron with the maximum number of synapses
    main_target_cell = target_synapse_counts.idxmax()
    total_synapses = target_synapse_counts.max()

    # To display the full equation, get the individual synapse counts from ASJL and ASJR
    # to the main target cell. We use .get(key, 0) to handle cases where a connection
    # from one of the pair might not exist.
    asjl_synapses_map = asj_connections[asj_connections['pre'] == 'ASJL'].set_index('post')['synapses']
    synapses_from_asjl = asjl_synapses_map.get(main_target_cell, 0)

    asjr_synapses_map = asj_connections[asj_connections['pre'] == 'ASJR'].set_index('post')['synapses']
    synapses_from_asjr = asjr_synapses_map.get(main_target_cell, 0)

    # Print the final analysis and the equation
    print(f"The main projection target of ASJ axons is the interneuron '{main_target_cell}'.")
    print("\nThe breakdown of the synapse count is as follows:")
    print(f"ASJL -> {main_target_cell}: {synapses_from_asjl} synapses")
    print(f"ASJR -> {main_target_cell}: {synapses_from_asjr} synapses")
    print("---------------------------------------")
    print(f"Total: {synapses_from_asjl} + {synapses_from_asjr} = {int(total_synapses)} synapses")


if __name__ == "__main__":
    find_asj_target()
    print("\n<<<AIY>>>")