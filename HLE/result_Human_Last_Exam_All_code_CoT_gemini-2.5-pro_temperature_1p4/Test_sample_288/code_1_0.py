import collections

try:
    from owconnectome import connectome
except ImportError:
    print("Please install the 'owconnectome' library first by running: pip install owconnectome")
    exit()

def find_main_target():
    """
    Finds the main projection target for a given source neuron pair in the C. elegans connectome.
    """
    # Load the connectome data
    c = connectome.Connectome()

    # The ASJ neuron pair
    source_neurons = ['ASJL', 'ASJR']
    
    # A dictionary to store synapse counts for each target neuron
    target_synapse_counts = collections.defaultdict(int)

    print(f"Analyzing outgoing chemical synapses from {', '.join(source_neurons)}...\n")

    # Iterate through the source neurons to get their connections
    for source_neuron in source_neurons:
        # Get all outgoing chemical synapses for the current neuron
        connections = c.connections(pre_cell=source_neuron, kind='chemical')
        for pre_cell, post_cell, num_synapses, conn_type in connections:
            # Add the synapse count to the target neuron's total
            target_synapse_counts[post_cell] += num_synapses
            
    if not target_synapse_counts:
        print("No outgoing connections found for ASJ neurons.")
        return

    # Find the neuron with the maximum number of synapses
    main_target_neuron = max(target_synapse_counts, key=target_synapse_counts.get)
    max_synapses = target_synapse_counts[main_target_neuron]

    # Print the full list of targets for context, sorted by synapse count
    print("Total synapse counts from ASJ (ASJL + ASJR) to target neurons:")
    sorted_targets = sorted(target_synapse_counts.items(), key=lambda item: item[1], reverse=True)
    for target, count in sorted_targets:
        print(f"ASJ -> {target}: {count} synapses")
    
    print("\n----------------------------------------")
    print(f"The main projection target of ASJ axons is {main_target_neuron} with a total of {max_synapses} synapses.")
    print("----------------------------------------")
    
    return main_target_neuron

if __name__ == "__main__":
    main_target = find_main_target()
