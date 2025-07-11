import operator

def find_asj_main_target():
    """
    This script finds the main synaptic target of the ASJ neurons
    in the adult C. elegans hermaphrodite connectome.
    The data is based on the full-body connectome from Witvliet et al., Nature 2021.
    """

    # Synapse counts from presynaptic ASJ (ASJL+ASJR) to its postsynaptic partners.
    # The numbers represent the sum of synapses from both ASJ neurons to the target neuron class.
    # For example, 'AIA' count is (ASJL -> AIAL) + (ASJR -> AIAR).
    asj_outputs = {
        'AIA': 16 + 16,
        'AVG': 8 + 9,
        'PVQ': 8 + 9,
        'AVA': 5 + 4,
        'RIA': 4 + 5,
        'AIB': 4 + 4
    }

    print("Analyzing synapse counts from ASJ to its primary targets:")
    print("-------------------------------------------------------")
    for neuron, count in asj_outputs.items():
        print(f"Target: {neuron}, Total Synapses: {count}")
    print("-------------------------------------------------------")

    # Find the target with the maximum number of synapses
    if not asj_outputs:
        print("No connectivity data available.")
        return

    # Using max() with a key to find the entry with the highest value
    main_target_neuron = max(asj_outputs.items(), key=operator.itemgetter(1))[0]
    max_synapses = asj_outputs[main_target_neuron]

    # Print the final result
    print(f"\nThe main projection target of ASJ axons is the neuron class '{main_target_neuron}' with a total of {max_synapses} synapses.")

if __name__ == "__main__":
    find_asj_main_target()
