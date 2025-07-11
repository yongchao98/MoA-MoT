def solve_mineral_ordering():
    """
    Identifies mineral specimens and orders them by increasing silicate structure complexity.
    """
    # Step 1 & 2: Define specimens with their identification and a complexity rank.
    # The rank is based on the silicate structure:
    # 1: Nesosilicate (Isolated Tetrahedra)
    # 2: Inosilicate (Chains)
    # 3: Phyllosilicate (Sheets)
    # 4: Tectosilicate (Framework)
    specimens = [
        {'label': 'A', 'mineral': 'Pyroxene', 'structure': 'Inosilicate (Single Chain)', 'complexity': 2},
        {'label': 'B', 'mineral': 'Mica', 'structure': 'Phyllosilicate (Sheet)', 'complexity': 3},
        {'label': 'C', 'mineral': 'Olivine', 'structure': 'Nesosilicate (Isolated Tetrahedra)', 'complexity': 1},
        {'label': 'D', 'mineral': 'Quartz', 'structure': 'Tectosilicate (Framework)', 'complexity': 4}
    ]

    print("--- Specimen Identification ---")
    # Create a dictionary for easy lookup
    specimen_dict = {s['label']: s for s in specimens}
    for label in sorted(specimen_dict.keys()):
        spec = specimen_dict[label]
        print(f"Specimen {spec['label']}: Identified as {spec['mineral']}, a {spec['structure']}.")

    # Step 3: Sort specimens by their complexity rank
    sorted_specimens = sorted(specimens, key=lambda x: x['complexity'])

    # Step 4: Extract the labels into the final ordered list
    ordered_labels = [s['label'] for s in sorted_specimens]

    print("\n--- Order of Increasing Silicate Structure Complexity ---")
    print(" -> ".join(ordered_labels))

solve_mineral_ordering()