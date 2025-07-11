def solve_silicate_complexity():
    """
    Identifies mineral specimens, classifies their silicate structures,
    and orders them by increasing structural complexity.
    """
    # Step 1: Define the mineral data based on visual identification and geological knowledge.
    # A complexity score is assigned: 1 (least complex) to 4 (most complex).
    specimens = [
        {
            "label": "A",
            "mineral": "Pyroxene",
            "structure_type": "Inosilicate (Single Chain)",
            "complexity_score": 2
        },
        {
            "label": "B",
            "mineral": "Mica",
            "structure_type": "Phyllosilicate (Sheet)",
            "complexity_score": 3
        },
        {
            "label": "C",
            "mineral": "Olivine",
            "structure_type": "Nesosilicate (Isolated Tetrahedra)",
            "complexity_score": 1
        },
        {
            "label": "D",
            "mineral": "Quartz",
            "structure_type": "Tectosilicate (Framework)",
            "complexity_score": 4
        }
    ]

    # Step 2: Sort the specimens based on their complexity score in increasing order.
    sorted_specimens = sorted(specimens, key=lambda s: s['complexity_score'])

    # Step 3: Print the identification for each specimen.
    print("Specimen Identification and Classification:")
    # We iterate through the original (unsorted) list to match the image labels A, B, C, D.
    unsorted_specimens = sorted(specimens, key=lambda s: s['label'])
    for specimen in unsorted_specimens:
        print(f"Specimen {specimen['label']}: Identified as {specimen['mineral']}, a {specimen['structure_type']} with complexity score {specimen['complexity_score']}.")

    # Step 4: Print the final ordered list.
    print("\nFinal order based on increasing silicate structure complexity:")
    
    # Construct the final equation-like string
    ordered_labels = [s['label'] for s in sorted_specimens]
    ordered_scores = [str(s['complexity_score']) for s in sorted_specimens]
    
    # Print the labels with arrows
    print("Labels: " + " -> ".join(ordered_labels))
    # Print the corresponding complexity scores
    print("Complexity: " + " < ".join(ordered_scores))


solve_silicate_complexity()