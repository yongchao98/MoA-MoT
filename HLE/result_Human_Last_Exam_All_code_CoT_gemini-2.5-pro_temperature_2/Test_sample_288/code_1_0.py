import operator

def find_asj_main_target():
    """
    This function identifies the main postsynaptic target of ASJ axons
    in the adult C. elegans hermaphrodite based on synapse number.
    The data is based on the connectome data from White et al., 1986 and subsequent reconstructions.
    """
    # Synapse counts from presynaptic neuron ASJ to its postsynaptic targets
    asj_axon_targets = {
        'AIY': 14,
        'AUA': 5,
        'PVQ': 5,
        'ASG': 4,
        'ASK': 4,
        'RIA': 3,
        'PVT': 2,
        'AVG': 1,
        'SAA': 1
    }

    # Find the target with the maximum number of synapses
    if not asj_axon_targets:
        print("No connectivity data available.")
        return

    # Using max() with a key to find the dictionary key with the largest value
    main_target_cell = max(asj_axon_targets, key=asj_axon_targets.get)
    max_synapses = asj_axon_targets[main_target_cell]

    print(f"The main projection target of ASJ axons is the {main_target_cell} neuron.")
    print(f"This is based on the total number of synapses received.\n")
    print("Synapse count equation:")
    print(f"{main_target_cell} = {max_synapses} synapses\n")

    print("Detailed breakdown of synapses from ASJ to its major targets:")
    # Sort the dictionary by synapse count in descending order for a clear ranking
    sorted_targets = sorted(asj_axon_targets.items(), key=operator.itemgetter(1), reverse=True)
    
    for cell, count in sorted_targets:
        print(f"ASJ -> {cell}: {count} synapses")

if __name__ == '__main__':
    find_asj_main_target()
