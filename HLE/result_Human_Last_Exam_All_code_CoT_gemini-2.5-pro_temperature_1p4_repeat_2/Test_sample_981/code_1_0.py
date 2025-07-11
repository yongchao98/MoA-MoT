def classify_fusion_experiments():
    """
    This function classifies fusion experiments based on the properties they use
    to twist their magnetic fields and prints the result.
    """
    # The properties are:
    # A: driving a toroidal current
    # B: elongating flux surfaces and making them rotate poloidally
    # C: making the magnetic axis non-planar

    classification = {
        "Uses only property B": [],
        "Uses both B and C": ["LHD", "Wendelstein 7-X", "NCSX"],
        "Uses A, B, and C": ["Tokamaks", "Reversed Field Pinches"]
    }

    print("Classification of fusion experiments based on magnetic field twisting properties:")
    print("-" * 70)

    for category, experiments in classification.items():
        print(f"{category}:")
        if not experiments:
            print("- None")
        else:
            for exp in experiments:
                print(f"- {exp}")
        # Add a newline for better readability between categories
        print()

# Execute the classification and print the results.
classify_fusion_experiments()
<<<Only B: None; B and C: LHD, Wendelstein 7-X, NCSX; A, B, and C: Tokamaks, Reversed Field Pinches>>>