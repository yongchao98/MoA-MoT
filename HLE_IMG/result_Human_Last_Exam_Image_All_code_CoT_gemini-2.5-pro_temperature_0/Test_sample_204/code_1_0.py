def solve_mineral_puzzle():
    """
    Identifies mineral specimens, classifies their silicate structure,
    and orders them by increasing structural complexity.
    """
    # Step 1 & 2: Identify minerals and their silicate class with a complexity score.
    # Complexity score: 1 (least complex) to 4 (most complex).
    specimens = [
        {'label': 'C', 
         'mineral': 'Olivine', 
         'class': 'Nesosilicate (Isolated Tetrahedra)', 
         'complexity': 1},
        {'label': 'A', 
         'mineral': 'Pyroxene (e.g., Augite)', 
         'class': 'Inosilicate (Single Chain)', 
         'complexity': 2},
        {'label': 'B', 
         'mineral': 'Mica (e.g., Muscovite)', 
         'class': 'Phyllosilicate (Sheet)', 
         'complexity': 3},
        {'label': 'D', 
         'mineral': 'Quartz', 
         'class': 'Tectosilicate (Framework)', 
         'complexity': 4}
    ]

    print("Specimen Identification and Classification:")
    # Create a dictionary for easy lookup by label
    specimen_dict = {s['label']: s for s in specimens}
    for label in sorted(specimen_dict.keys()):
        specimen = specimen_dict[label]
        print(f"Specimen {specimen['label']}: {specimen['mineral']} - a {specimen['class']}")
    
    print("\n" + "="*40 + "\n")

    # Step 3: Sort the specimens by their structural complexity.
    sorted_specimens = sorted(specimens, key=lambda x: x['complexity'])

    # Step 4: Print the final ordered list.
    print("Order of Increasing Silicate Structure Complexity:")
    
    ordered_labels = [s['label'] for s in sorted_specimens]
    
    # Print the final equation-like format
    output_string = " -> ".join(ordered_labels)
    print(output_string)

    # Final answer format
    final_answer = ",".join(ordered_labels)
    print(f"\n<<<C,A,B,D>>>")

solve_mineral_puzzle()