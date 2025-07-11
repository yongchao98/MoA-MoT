def solve_mineral_complexity():
    """
    Identifies mineral specimens, classifies their silicate structure,
    and orders them by increasing structural complexity.
    """

    # Data for each specimen: Mineral Name, Silicate Structure, Complexity Rank
    specimens = {
        'A': {'name': 'Pyroxene', 'structure': 'Single-Chain Silicate (Inosilicate)', 'complexity': 2},
        'B': {'name': 'Mica/Chlorite', 'structure': 'Sheet Silicate (Phyllosilicate)', 'complexity': 3},
        'C': {'name': 'Olivine', 'structure': 'Isolated Tetrahedra (Nesosilicate)', 'complexity': 1},
        'D': {'name': 'Quartz', 'structure': 'Framework Silicate (Tectosilicate)', 'complexity': 4}
    }

    print("Step 1: Identify each specimen and its silicate structure.")
    for label, data in specimens.items():
        print(f"Specimen {label}: Identified as {data['name']}, which is a {data['structure']}.")

    print("\nStep 2: Order the specimens by increasing silicate structure complexity.")
    print("The order of complexity is: Isolated < Single-Chain < Sheet < Framework.\n")

    # Sort the specimens based on the complexity rank
    sorted_specimens = sorted(specimens.items(), key=lambda item: item[1]['complexity'])

    # Print the final ordered list
    print("Final Order (from least to most complex):")
    
    ordered_labels = []
    for label, data in sorted_specimens:
        ordered_labels.append(label)
        print(f"{label} ({data['name']}) - {data['structure']}")

    final_answer = ', '.join(ordered_labels)
    print(f"\nThe final order is: {final_answer}")
    
    # The final answer format required by the platform
    print(f"\n<<<{final_answer}>>>")

solve_mineral_complexity()