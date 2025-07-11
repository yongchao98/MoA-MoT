def classify_fusion_experiments():
    """
    Classifies fusion experiments based on their methods for twisting magnetic fields.

    Properties:
    (A) Driving a toroidal current
    (B) Elongating flux surfaces and poloidal rotation
    (C) Non-planar magnetic axis
    """

    # Data structure: Experiment -> Set of properties used
    experiments = {
        "Tokamaks": {"A", "B"},
        "LHD": {"B", "C"},
        "Wendelstein 7-X": {"B", "C"},
        "NCSX": {"A", "B", "C"},
        "Reversed Field Pinches": {"A", "B"}
    }

    # Groups to hold the categorized experiments
    only_b = []
    b_and_c = []
    a_b_and_c = []
    # An additional category is needed for a complete classification
    a_and_b = []

    # Iterate and classify each experiment
    for name, properties in experiments.items():
        if properties == {"B"}:
            only_b.append(name)
        elif properties == {"B", "C"}:
            b_and_c.append(name)
        elif properties == {"A", "B", "C"}:
            a_b_and_c.append(name)
        elif properties == {"A", "B"}:
            a_and_b.append(name)

    # Print the final results in a structured format
    print("Classification of Fusion Experiments:")
    print("-" * 40)

    print("Use only property B:")
    if not only_b:
        print(" - None")
    else:
        for exp in only_b:
            print(f" - {exp}")

    print("\nUse both properties B and C:")
    if not b_and_c:
        print(" - None")
    else:
        for exp in b_and_c:
            print(f" - {exp}")

    print("\nUse properties A, B, and C:")
    if not a_b_and_c:
        print(" - None")
    else:
        for exp in a_b_and_c:
            print(f" - {exp}")

    # Although not explicitly requested, this category is required to place all listed experiments
    print("\nUse properties A and B:")
    if not a_and_b:
        print(" - None")
    else:
        for exp in a_and_b:
            print(f" - {exp}")

classify_fusion_experiments()