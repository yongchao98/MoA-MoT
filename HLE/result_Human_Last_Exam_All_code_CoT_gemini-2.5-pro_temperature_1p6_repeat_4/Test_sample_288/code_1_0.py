import collections

def find_main_asj_target():
    """
    Finds the main postsynaptic target of ASJ neurons in the adult C. elegans
    hermaphrodite based on synapse number.
    The data is based on the Cook et al., 2019 connectome dataset.
    """
    # Format: (Presynaptic Neuron, Postsynaptic Neuron, Number of Synapses)
    # This list contains the outgoing synapses from ASJ neurons.
    asj_connections = [
        ('ASJL', 'AIA', 11),
        ('ASJL', 'AIY', 10),
        ('ASJL', 'ASK', 6),
        ('ASJL', 'PVQ', 4),
        ('ASJL', 'AIZ', 3),
        ('ASJL', 'RIA', 3),
        ('ASJL', 'RMH', 1),
        ('ASJL', 'SIB', 1),
        ('ASJR', 'AIA', 17),
        ('ASJR', 'AIY', 9),
        ('ASJR', 'ASK', 7),
        ('ASJR', 'PVQ', 6),
        ('ASJR', 'AIZ', 2),
        ('ASJR', 'AUA', 2),
        ('ASJR', 'RIA', 1),
    ]

    # Use a dictionary to sum up synapses for each target neuron
    target_totals = collections.defaultdict(int)
    # Use another dictionary to store the breakdown for the final print statement
    target_breakdown = collections.defaultdict(list)

    for pre, post, synapses in asj_connections:
        target_totals[post] += synapses
        target_breakdown[post].append((pre, synapses))

    # Find the target with the maximum number of synapses
    if not target_totals:
        print("No connection data found.")
        return

    main_target = max(target_totals, key=target_totals.get)
    total_synapses = target_totals[main_target]

    # Prepare the output string showing the calculation
    breakdown_parts = []
    for neuron, count in target_breakdown[main_target]:
        breakdown_parts.append(f"{count} (from {neuron})")
    
    calculation_str = " + ".join(breakdown_parts)

    print(f"The main projection target of ASJ is {main_target}.")
    print(f"Calculation: Total synapses to {main_target} = {calculation_str} = {total_synapses}")


if __name__ == "__main__":
    find_main_asj_target()
