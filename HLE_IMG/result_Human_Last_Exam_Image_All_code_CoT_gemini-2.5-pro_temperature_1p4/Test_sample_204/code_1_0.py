def solve_silicate_complexity():
    """
    Identifies mineral specimens and orders them by increasing silicate structure complexity.
    """
    # Step 1: Identify each specimen and its silicate class.
    specimens = {
        'A': {'mineral': 'Pyroxene', 'class': 'Inosilicate (Single Chain)', 'complexity': 2},
        'B': {'mineral': 'Mica/Chlorite', 'class': 'Phyllosilicate (Sheet)', 'complexity': 3},
        'C': {'mineral': 'Olivine', 'class': 'Nesosilicate (Isolated Tetrahedra)', 'complexity': 1},
        'D': {'mineral': 'Quartz', 'class': 'Tectosilicate (Framework)', 'complexity': 4}
    }

    print("Specimen Identification and Silicate Class:")
    for letter, data in sorted(specimens.items()):
        print(f"Specimen {letter}: {data['mineral']} ({data['class']})")

    # Step 2: Sort the specimens based on their complexity index.
    # The complexity index is pre-assigned: 1 (neso), 2 (ino), 3 (phyllo), 4 (tecto).
    sorted_specimens = sorted(specimens.items(), key=lambda item: item[1]['complexity'])
    
    # Step 3: Extract the letters for the final ordered list.
    ordered_letters = [item[0] for item in sorted_specimens]

    print("\nOrder of increasing silicate structure complexity:")
    # Using ' -> ' to visually represent the order.
    print(" -> ".join(ordered_letters))

solve_silicate_complexity()
<<<C -> A -> B -> D>>>