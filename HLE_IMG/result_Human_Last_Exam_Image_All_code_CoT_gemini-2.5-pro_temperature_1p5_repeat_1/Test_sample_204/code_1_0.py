def solve_mineral_puzzle():
    """
    Identifies mineral specimens, classifies their silicate structures,
    and orders them by increasing structural complexity.
    """
    # Step 1-4: Identify specimens and classify their silicate structures.
    # The complexity is ranked based on the degree of silica tetrahedra polymerization:
    # 1: Nesosilicate (isolated tetrahedra, SiO4)
    # 2: Inosilicate (single chains, SiO3)
    # 3: Phyllosilicate (sheets, Si2O5)
    # 4: Tectosilicate (framework, SiO2)
    specimens = {
        'C': {'mineral': 'Olivine', 'structure': 'Nesosilicate (Isolated Tetrahedra)', 'complexity': 1},
        'A': {'mineral': 'Pyroxene', 'structure': 'Inosilicate (Single Chain)', 'complexity': 2},
        'B': {'mineral': 'Mica', 'structure': 'Phyllosilicate (Sheet)', 'complexity': 3},
        'D': {'mineral': 'Quartz', 'structure': 'Tectosilicate (Framework)', 'complexity': 4}
    }

    print("Specimen Identification and Silicate Structure:\n")
    # Print details for each specimen in alphabetical order first.
    for letter in sorted(specimens.keys()):
        info = specimens[letter]
        print(f"Specimen {letter}:")
        print(f"  - Identification: {info['mineral']}")
        print(f"  - Silicate Structure: {info['structure']}")
        print(f"  - Complexity Rank: {info['complexity']}\n")

    # Step 5: Sort the specimens by their structural complexity.
    sorted_specimens = sorted(specimens.items(), key=lambda item: item[1]['complexity'])

    # Extract the letters for the final ordered list
    ordered_list = [item[0] for item in sorted_specimens]

    print("--------------------------------------------------")
    print("Order of Increasing Silicate Structure Complexity:")
    print(" -> ".join(ordered_list))

solve_mineral_puzzle()