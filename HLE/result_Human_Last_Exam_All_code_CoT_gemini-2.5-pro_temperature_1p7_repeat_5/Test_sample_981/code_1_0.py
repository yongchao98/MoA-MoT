def classify_fusion_experiments():
    """
    Classifies fusion experiments based on their method of twisting magnetic fields.

    The properties are:
    (A) Driving a toroidal current.
    (B) Elongating the flux surfaces and making them rotate poloidally.
    (C) Making the magnetic axis non-planar.

    The classification logic is as follows:
    - "Only property B": Planar-axis devices where the twist is fundamental.
       This category includes Tokamaks and RFPs, where property (A) is the mechanism
       to create property (B), and property (C) is absent.
    - "Both B and C": Pure stellarators with no significant current (A), which rely on
       3D coil shaping (creating B) and have a non-planar axis (C).
    - "A, B, and C": Hybrid devices that use both a toroidal current (A) and
       3D shaping with a non-planar axis (B and C).
    """

    # Mapping each experiment to its classification
    experiments_classification = {
        "Tokamaks": "Only property B",
        "Reversed Field Pinches": "Only property B",
        "LHD": "Both B and C",
        "Wendelstein 7-X": "Both B and C",
        "NCSX": "A, B, and C"
    }

    # Group experiments by category for organized printing
    grouped_results = {
        "Only property B": [],
        "Both B and C": [],
        "A, B, and C": []
    }

    for exp, category in experiments_classification.items():
        if category in grouped_results:
            grouped_results[category].append(exp)

    # Print the results in a clear format
    print("Classification of Fusion Experiments:")
    print("-" * 40)

    # Print experiments using only property B
    category_b = "Only property B"
    print(f"{category_b}:")
    if grouped_results[category_b]:
        for experiment in grouped_results[category_b]:
            print(f"  - {experiment}")
    else:
        print("  - None")
    print()

    # Print experiments using both B and C
    category_bc = "Both B and C"
    print(f"{category_bc}:")
    if grouped_results[category_bc]:
        for experiment in grouped_results[category_bc]:
            print(f"  - {experiment}")
    else:
        print("  - None")
    print()
    
    # Print experiments using A, B, and C
    category_abc = "A, B, and C"
    print(f"{category_abc}:")
    if grouped_results[category_abc]:
        for experiment in grouped_results[category_abc]:
            print(f"  - {experiment}")
    else:
        print("  - None")

# Execute the classification function
classify_fusion_experiments()