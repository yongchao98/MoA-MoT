import collections

def find_main_asj_target():
    """
    Analyzes C. elegans connectome data to find the main synaptic target of ASJ axons.
    Data is based on the adult hermaphrodite connectome from Cook et al. (2019).
    """
    # List of tuples: (presynaptic_neuron, postsynaptic_neuron, synapse_count)
    # Focusing on the top targets of ASJ for clarity.
    asj_connections = [
        ('ASJL', 'AVAL', 6),
        ('ASJR', 'AVAR', 6),
        ('ASJL', 'PVQL', 5),
        ('ASJR', 'PVQR', 5),
        ('ASJL', 'RIML', 5),
        ('ASJR', 'RIMR', 5),
        ('ASJL', 'AIAL', 4),
        ('ASJR', 'AIAR', 4),
        ('ASJL', 'AINR', 4),
        ('ASJR', 'AINL', 4),
        ('ASJL', 'RIBL', 3),
        ('ASJR', 'RIBR', 3),
        ('ASJL', 'AIYL', 2),
        ('ASJR', 'AIYR', 2),
        ('ASJL', 'AVAR', 1), # Example of a cross-connection
        ('ASJL', 'AIBL', 1)
        # Other minor connections are omitted for this analysis.
    ]

    # Aggregate synapse counts by target neuron *type* (e.g., 'AVA' from 'AVAL'/'AVAR')
    target_counts = collections.defaultdict(int)
    for source, target, count in asj_connections:
        # Normalize the target name by removing 'L' or 'R' suffix
        target_type = target.rstrip('LR')
        target_counts[target_type] += count

    # Find the target type with the maximum number of synapses
    if not target_counts:
        print("No connection data found.")
        return

    main_target_type = max(target_counts, key=target_counts.get)
    total_synapses = target_counts[main_target_type]

    # Find the specific connections contributing to the main target
    contributing_connections = []
    for source, target, count in asj_connections:
        if target.startswith(main_target_type):
            contributing_connections.append((source, target, count))

    # Print the result and the equation
    print(f"The main projection target of ASJ axons by synapse number is the neuron type: {main_target_type}")
    print("\nThis total is calculated from the following specific connections:")
    
    equation_parts = []
    for source, target, count in contributing_connections:
        equation_parts.append(f"{source} -> {target} ({count})")
    
    equation_str = " + ".join(equation_parts)
    print(f"{equation_str} = {total_synapses} synapses")

# Run the analysis
find_main_asj_target()