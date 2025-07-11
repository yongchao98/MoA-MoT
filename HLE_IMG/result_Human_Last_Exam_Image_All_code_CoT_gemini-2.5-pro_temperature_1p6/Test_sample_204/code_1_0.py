def solve_mineral_ordering():
    """
    Identifies mineral specimens, classifies their silicate structure,
    and orders them by increasing structural complexity.
    """
    # Step 1: Define the specimens with their identified properties.
    # Complexity is ranked from 1 (least complex) to 4 (most complex).
    specimens = [
        {
            "label": "A",
            "mineral": "Pyroxene/Amphibole",
            "structure": "Inosilicate (Chain)",
            "complexity_rank": 2
        },
        {
            "label": "B",
            "mineral": "Mica/Chlorite",
            "structure": "Phyllosilicate (Sheet)",
            "complexity_rank": 3
        },
        {
            "label": "C",
            "mineral": "Olivine",
            "structure": "Nesosilicate (Isolated Tetrahedra)",
            "complexity_rank": 1
        },
        {
            "label": "D",
            "mineral": "Quartz",
            "structure": "Tectosilicate (Framework)",
            "complexity_rank": 4
        }
    ]

    # Step 2: Sort the specimens based on their complexity rank.
    sorted_specimens = sorted(specimens, key=lambda x: x['complexity_rank'])

    # Step 3: Print the identification and ordering logic.
    print("Identification and Ordering by Silicate Structure Complexity:")
    print("-" * 60)
    for specimen in sorted_specimens:
        print(f"Specimen {specimen['label']}:")
        print(f"  - Identification: {specimen['mineral']}")
        print(f"  - Structure Type: {specimen['structure']}")
        print(f"  - Complexity Rank: {specimen['complexity_rank']}")
    print("-" * 60)

    # Step 4: Construct and print the final answer string.
    final_order = ", ".join([s['label'] for s in sorted_specimens])
    print(f"The specimens ordered by increasing silicate structure complexity are: {final_order}")

solve_mineral_ordering()
<<<C, A, B, D>>>