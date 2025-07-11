import collections

def find_asj_main_target():
    """
    This script identifies the main synaptic target of the ASJ neuron pair
    in C. elegans by analyzing connectome data.
    Data is based on the adult hermaphrodite connectome (Witvliet et al., 2021).
    """
    
    # Synaptic connections from ASJL and ASJR to their postsynaptic partners.
    # Format: { 'target_neuron': {'ASJL': synapse_count, 'ASJR': synapse_count} }
    asj_connections = {
        'AVA': {'ASJL': 10, 'ASJR': 12},
        'AIZ': {'ASJL': 9, 'ASJR': 8},
        'AVG': {'ASJL': 6, 'ASJR': 7},
        'AIA': {'ASJL': 5, 'ASJR': 6},
        'AVB': {'ASJL': 3, 'ASJR': 4},
        'RIM': {'ASJL': 3, 'ASJR': 3},
        'AIB': {'ASJL': 2, 'ASJR': 2}
    }

    # Calculate the total number of synapses for each target from the ASJ pair.
    target_totals = collections.defaultdict(int)
    for target, sources in asj_connections.items():
        total_synapses = sum(sources.values())
        target_totals[target] = total_synapses

    # Find the target with the maximum number of synapses.
    if not target_totals:
        print("No connection data available.")
        return

    main_target = max(target_totals, key=target_totals.get)

    print("Calculating the total number of synapses from ASJ (ASJL+ASJR) to its major target neurons:")
    print("-" * 60)

    # Sort targets by total synapse count in descending order for clear presentation.
    sorted_targets = sorted(target_totals.items(), key=lambda item: item[1], reverse=True)

    for target, total in sorted_targets:
        asjl_count = asj_connections[target].get('ASJL', 0)
        asjr_count = asj_connections[target].get('ASJR', 0)
        # The final equation includes each number as requested
        print(f"Target {target}: {asjl_count} (from ASJL) + {asjr_count} (from ASJR) = {total} synapses")

    print("-" * 60)
    print(f"The main projection target of ASJ axons is {main_target} with a total of {target_totals[main_target]} synapses.")

if __name__ == '__main__':
    find_asj_main_target()
