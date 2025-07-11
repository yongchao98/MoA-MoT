def solve_mineral_complexity():
    """
    Identifies mineral specimens, classifies their silicate structures,
    and orders them by increasing structural complexity.
    """
    # Step 1: Define the specimens with their identified mineral,
    # silicate structure class, and a complexity rank.
    # Complexity rank: 1 (least complex) to 4 (most complex)
    specimens = [
        {'id': 'A', 'mineral': 'Pyroxene', 'structure': 'Single-chain silicate (Inosilicate)', 'complexity': 2},
        {'id': 'B', 'mineral': 'Mica/Chlorite', 'structure': 'Sheet silicate (Phyllosilicate)', 'complexity': 3},
        {'id': 'C', 'mineral': 'Olivine', 'structure': 'Isolated tetrahedra (Nesosilicate)', 'complexity': 1},
        {'id': 'D', 'mineral': 'Quartz', 'structure': 'Framework silicate (Tectosilicate)', 'complexity': 4}
    ]

    print("--- Specimen Identification and Classification ---")
    for specimen in specimens:
        print(f"Specimen {specimen['id']}:")
        print(f"  - Identification: {specimen['mineral']}")
        print(f"  - Silicate Structure: {specimen['structure']}")
        print("-" * 20)

    # Step 2: Sort the specimens based on their complexity rank
    sorted_specimens = sorted(specimens, key=lambda x: x['complexity'])

    # Step 3: Generate the final ordered list of specimen IDs
    ordered_ids = [s['id'] for s in sorted_specimens]

    print("\n--- Final Order by Increasing Silicate Structure Complexity ---")
    # The final print statement is formatted to show the sequence clearly.
    print(f"The order is: {ordered_ids[0]}, {ordered_ids[1]}, {ordered_ids[2]}, {ordered_ids[3]}")

# Run the function to get the solution
solve_mineral_complexity()
