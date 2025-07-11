def classify_fusion_experiments():
    """
    Classifies fusion experiments based on their magnetic field twisting properties.
    Properties:
    A: Driving a toroidal current
    B: Elongating flux surfaces (shaping)
    C: Non-planar magnetic axis
    """
    # Database of experiments and their properties
    experiments = {
        "Tokamaks": {'A', 'B'},
        "LHD": {'B', 'C'},
        "Wendelstein 7-X": {'B', 'C'},
        "NCSX": {'A', 'B', 'C'},
        "Reversed Field Pinches": {'A', 'B'}
    }

    # The categories requested by the user
    categories_to_find = {
        "Only B": {'B'},
        "B and C": {'B', 'C'},
        "A, B, and C": {'A', 'B', 'C'}
    }

    # Store the results to be printed
    results = []

    # Classify experiments into the requested categories
    for category_name, properties_set in categories_to_find.items():
        # Find all experiments whose set of properties exactly matches the category's set
        matching_experiments = sorted([
            name for name, props in experiments.items() if props == properties_set
        ])

        # Format the result string for the current category
        if not matching_experiments:
            result_line = f"{category_name}: None"
        else:
            result_line = f"{category_name}: {', '.join(matching_experiments)}"
        results.append(result_line)

    # Print the final, formatted classification
    final_output = "\n".join(results)
    print(final_output)

    # Also provide the answer in the special format for direct extraction
    print(f"\n<<<Only B: None\nB and C: LHD, Wendelstein 7-X\nA, B, and C: NCSX>>>")

# Execute the function
classify_fusion_experiments()