def solve_mineral_complexity():
    """
    Identifies mineral specimens and orders them by increasing silicate structure complexity.
    """
    # Step 1 & 2: Identify minerals and their silicate structures
    # The complexity is ranked from 1 (least complex) to 4 (most complex).
    minerals = [
        {'label': 'C', 'name': 'Olivine', 'structure': 'Nesosilicate (Isolated Tetrahedra)', 'complexity': 1},
        {'label': 'A', 'name': 'Pyroxene', 'structure': 'Inosilicate (Single Chain)', 'complexity': 2},
        {'label': 'B', 'name': 'Mica/Chlorite', 'structure': 'Phyllosilicate (Sheet)', 'complexity': 3},
        {'label': 'D', 'name': 'Quartz', 'structure': 'Tectosilicate (Framework)', 'complexity': 4}
    ]

    # Step 3: Sort the minerals based on their complexity rank
    sorted_minerals = sorted(minerals, key=lambda x: x['complexity'])

    # Step 4: Print the identification for each specimen
    print("Specimen Identification and Silicate Structure:")
    # Print the details in the original A, B, C, D order for clarity
    for mineral in sorted(minerals, key=lambda x: x['label']):
        print(f"  Specimen {mineral['label']}: {mineral['name']} - {mineral['structure']}")

    # Step 5: Print the final ordered list
    print("\nOrder of increasing silicate structure complexity:")
    
    ordered_labels = [mineral['label'] for mineral in sorted_minerals]
    
    # Print the equation-like format
    print(f"{ordered_labels[0]} < {ordered_labels[1]} < {ordered_labels[2]} < {ordered_labels[3]}")

solve_mineral_complexity()