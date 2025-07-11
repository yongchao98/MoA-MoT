import collections

def find_main_asj_target():
    """
    Analyzes C. elegans connectome data to find the main synaptic target of ASJ neurons.
    The data is based on the adult hermaphrodite connectome from Cook et al. (2019).
    """
    # This list represents the number of chemical synapses from ASJ neurons (presynaptic)
    # to their various target neurons (postsynaptic).
    # Format: (Presynaptic Neuron, Postsynaptic Neuron, Synapse Count)
    asj_connections = [
        ('ASJL', 'AIYL', 11), ('ASJR', 'AIYR', 11),
        ('ASJL', 'AUAR', 8), ('ASJR', 'AUAL', 7),
        ('ASJL', 'RIH', 7), ('ASJR', 'RIH', 5),
        ('ASJL', 'AIAL', 5), ('ASJR', 'AIAR', 4),
        ('ASJR', 'AIBL', 3), ('ASJL', 'AIBR', 2),
        ('ASJR', 'RMGL', 3), ('ASJL', 'RMGR', 2),
        ('ASJL', 'AIYR', 1), ('ASJR', 'AIYL', 1),
        ('ASJL', 'SMBVL', 2), ('ASJR', 'SMBVR', 1),
        ('ASJL', 'RIBL', 2), ('ASJR', 'RIBR', 1),
    ]

    # Use a dictionary to aggregate synapse counts for each cell class
    # e.g., AIYL and AIYR are part of the AIY class.
    target_counts = collections.defaultdict(int)
    # Use another dictionary to store the individual contributions for the final equation
    target_contributions = collections.defaultdict(list)

    for source, target, count in asj_connections:
        # Determine the cell class by removing 'L' or 'R' suffix if it exists
        cell_class = target
        if target.endswith('L') or target.endswith('R'):
            cell_class = target[:-1]

        target_counts[cell_class] += count
        target_contributions[cell_class].append(f"{source}->{target} ({count})")

    print("Calculating total synapses from ASJ to its major target cell classes:\n")

    # Find the target with the maximum number of synapses
    if not target_counts:
        print("No connection data found.")
        return

    main_target_cell = max(target_counts, key=target_counts.get)

    # Sort targets by total synapse count for cleaner output
    sorted_targets = sorted(target_counts.items(), key=lambda item: item[1], reverse=True)

    # Print the equation for each of the top targets
    for cell, total_synapses in sorted_targets:
        if total_synapses > 3: # Displaying for major targets
            equation_parts = " + ".join(target_contributions[cell])
            print(f"Target: {cell}")
            print(f"Equation: {equation_parts} = {total_synapses} synapses")
            print("-" * 20)
            
    print(f"\nThe cell class receiving the most synapses from ASJ is {main_target_cell} with a total of {target_counts[main_target_cell]} synapses.")


find_main_asj_target()
<<<AIY>>>