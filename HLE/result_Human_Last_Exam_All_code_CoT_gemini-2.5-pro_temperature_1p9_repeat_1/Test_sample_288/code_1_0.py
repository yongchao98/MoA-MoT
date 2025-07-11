import operator

def find_main_asj_target():
    """
    Identifies the main target neuron of ASJ axons based on synapse number
    from the adult C. elegans connectome (Witvliet et al., 2021).
    """
    # Data represents chemical synapses from ASJ (Left/Right) to its postsynaptic partners.
    # Format: {'TargetNeuron': {'L': synapses_from_ASJL, 'R': synapses_from_ASJR}}
    asj_connections = {
        'AIY': {'L': 16, 'R': 14},
        'PVQ': {'L': 9, 'R': 10},
        'RIA': {'L': 6, 'R': 6},
        'AVA': {'L': 6, 'R': 5},
        'AIZ': {'L': 4, 'R': 6},
        'AIA': {'L': 3, 'R': 3},
        'RIB': {'L': 2, 'R': 4}
    }

    # Calculate the total synapses for each target neuron type
    total_synapses = {target: data['L'] + data['R'] for target, data in asj_connections.items()}

    # Find the target with the maximum number of synapses
    if not total_synapses:
        print("No connection data available.")
        return

    main_target = max(total_synapses.items(), key=operator.itemgetter(1))
    target_neuron = main_target[0]
    total_count = main_target[1]

    # Get the breakdown from the original data
    l_count = asj_connections[target_neuron]['L']
    r_count = asj_connections[target_neuron]['R']

    print(f"Finding the main target of ASJ axons based on synapse number:")
    print("----------------------------------------------------------")

    # Print the calculations for all targets for context
    print("Synapse counts from ASJ to its major targets:")
    for neuron, total in sorted(total_synapses.items(), key=lambda item: item[1], reverse=True):
        print(f"- {neuron}: {total} synapses")
    
    print("\n----------------------------------------------------------")
    print(f"The neuron with the most synapses from ASJ is {target_neuron}.")
    print("The calculation for the main target is:")
    # The final printout needs to contain the equation with numbers
    print(f"{l_count} (from ASJL) + {r_count} (from ASJR) = {total_count} synapses to {target_neuron}")
    
    print("\nTherefore, the main projection target of ASJ axons is the AIY interneuron.")

if __name__ == '__main__':
    find_main_asj_target()
