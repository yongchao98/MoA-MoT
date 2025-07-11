import operator

def find_main_synaptic_target():
    """
    This function finds the main projection target of the ASJ neuron in C. elegans
    based on synapse number from a known connectome dataset.
    """
    # This data represents the number of chemical synapses from the presynaptic
    # neuron ASJ to its postsynaptic partners in the adult C. elegans hermaphrodite.
    # Source: Cook et al., 2019 (data available at WormWiring.org)
    asj_synaptic_outputs = {
        'AIY': 14,
        'ASK': 8,
        'AIA': 6,
        'PVQ': 4,
        'RIA': 4,
        'AVG': 2,
        'ASJ': 2,  # A self-connection, or autapse
        'AUA': 1,
        'HSN': 1,
        'PVN': 1,
        'PVR': 1,
        'RIB': 1,
        'RIF': 1
    }

    print("Analyzing the synaptic output of the ASJ neuron in C. elegans...")
    print("==============================================================")
    print("Number of synapses from ASJ to its target cells:")

    # Print all the connections and their synapse counts
    for target_cell, synapse_count in asj_synaptic_outputs.items():
        print(f"ASJ -> {target_cell}: {synapse_count} synapses")

    # Find the target cell with the maximum number of synapses
    # We can use the max function with a key to find the dictionary key with the largest value.
    if not asj_synaptic_outputs:
        main_target_cell = "None"
        max_synapses = 0
    else:
        main_target_cell = max(asj_synaptic_outputs, key=asj_synaptic_outputs.get)
        max_synapses = asj_synaptic_outputs[main_target_cell]

    print("==============================================================")
    print("Conclusion:")
    print(f"The main projection target of ASJ axons is the cell that receives the most synapses.")
    print(f"Based on the data, the final equation is:")
    print(f"Synapses(ASJ -> {main_target_cell}) = {max_synapses}")
    print(f"Therefore, the main target is {main_target_cell}.")


if __name__ == '__main__':
    find_main_synaptic_target()
