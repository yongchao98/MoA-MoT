def solve_silicate_puzzle():
    """
    Identifies mineral specimens and orders them by increasing silicate structure complexity.
    """
    specimens = [
        {
            "label": "A",
            "mineral": "Pyroxene (e.g., Augite)",
            "structure": "Inosilicate (Single Chain)",
            "complexity": 2
        },
        {
            "label": "B",
            "mineral": "Mica / Chlorite",
            "structure": "Phyllosilicate (Sheet)",
            "complexity": 3
        },
        {
            "label": "C",
            "mineral": "Olivine",
            "structure": "Nesosilicate (Isolated Tetrahedra)",
            "complexity": 1
        },
        {
            "label": "D",
            "mineral": "Quartz",
            "structure": "Tectosilicate (Framework)",
            "complexity": 4
        }
    ]

    print("Specimen Identification:")
    # Print details for each specimen in alphabetical order
    for specimen in sorted(specimens, key=lambda x: x['label']):
        print(f"Specimen {specimen['label']}: {specimen['mineral']}, a {specimen['structure']}.")

    # Sort the specimens by complexity
    sorted_specimens = sorted(specimens, key=lambda x: x['complexity'])

    # Get the ordered labels
    ordered_labels = [s['label'] for s in sorted_specimens]

    print("\nOrder of increasing silicate structure complexity:")
    # The 'end' parameter is used to print the elements on the same line separated by ' -> '
    print(" -> ".join(ordered_labels))

solve_silicate_puzzle()