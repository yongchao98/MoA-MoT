def solve_mineral_complexity():
    """
    Identifies mineral specimens and orders them by increasing silicate structure complexity.
    """
    # Step 1 & 2: Identify specimens and classify their silicate structure.
    # We assign a complexity score: 1 (simplest) to 4 (most complex).
    specimens = [
        {'label': 'C', 'mineral': 'Olivine', 'group': 'Nesosilicate (Isolated Tetrahedra)', 'complexity': 1},
        {'label': 'A', 'mineral': 'Pyroxene', 'group': 'Inosilicate (Single Chain)', 'complexity': 2},
        {'label': 'B', 'mineral': 'Mica/Chlorite', 'group': 'Phyllosilicate (Sheet)', 'complexity': 3},
        {'label': 'D', 'mineral': 'Quartz', 'group': 'Tectosilicate (Framework)', 'complexity': 4}
    ]

    # Step 3: Sort the specimens based on the complexity score.
    sorted_specimens = sorted(specimens, key=lambda x: x['complexity'])

    # Step 4: Print the identification and the final ordered list.
    print("Identifying and ordering mineral specimens by silicate structure complexity:\n")
    for specimen in sorted_specimens:
        print(f"Specimen {specimen['label']}:")
        print(f"  - Identification: {specimen['mineral']}")
        print(f"  - Silicate Group: {specimen['group']}")
        print(f"  - Complexity Level: {specimen['complexity']}\n")

    ordered_labels = [s['label'] for s in sorted_specimens]
    final_order_str = ", ".join(ordered_labels)

    print(f"The order of specimens in increasing silicate structure complexity is: {final_order_str}")

solve_mineral_complexity()
<<<C, A, B, D>>>