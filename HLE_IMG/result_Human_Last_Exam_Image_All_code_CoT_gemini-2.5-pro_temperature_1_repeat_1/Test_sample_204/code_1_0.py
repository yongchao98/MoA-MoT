def solve_mineral_complexity():
    """
    Identifies mineral specimens, classifies their silicate structure,
    and orders them by increasing structural complexity.
    """

    # Data for each specimen based on visual identification and geological knowledge.
    # Complexity is ranked from 1 (least complex) to 4 (most complex).
    specimens = [
        {'label': 'A', 'mineral': 'Pyroxene', 'structure': 'Inosilicate (Single Chain)', 'complexity_rank': 2},
        {'label': 'B', 'mineral': 'Mica', 'structure': 'Phyllosilicate (Sheet)', 'complexity_rank': 3},
        {'label': 'C', 'mineral': 'Olivine', 'structure': 'Nesosilicate (Isolated Tetrahedra)', 'complexity_rank': 1},
        {'label': 'D', 'mineral': 'Quartz', 'structure': 'Tectosilicate (Framework)', 'complexity_rank': 4}
    ]

    # Sort the specimens based on their complexity rank
    sorted_specimens = sorted(specimens, key=lambda x: x['complexity_rank'])

    print("Identification of Specimens and their Silicate Structures:")
    print("-" * 60)
    for specimen in specimens:
        print(f"Specimen {specimen['label']}: Identified as {specimen['mineral']}, a {specimen['structure']}.")
    print("-" * 60)

    # Get the final ordered list of labels
    ordered_labels = [s['label'] for s in sorted_specimens]

    print("\nSpecimens ordered by increasing silicate structure complexity:")
    # The prompt asks to "output each number in the final equation", which is interpreted here
    # as printing the final sequence clearly.
    final_order_string = ", ".join(ordered_labels)
    print(final_order_string)

solve_mineral_complexity()
<<<C, A, B, D>>>