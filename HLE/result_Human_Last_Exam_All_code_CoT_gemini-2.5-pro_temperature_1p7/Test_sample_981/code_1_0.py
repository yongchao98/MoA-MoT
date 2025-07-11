def classify_fusion_experiments():
    """
    Classifies fusion experiments based on their use of three magnetic field twisting properties.

    Property A: Driving a toroidal current.
    Property B: Elongating flux surfaces (shaping/external transform).
    Property C: Making the magnetic axis non-planar.

    The classification logic is derived from the standard understanding of these devices,
    mapped onto the specific categories provided in the problem.
    """

    # Data represents the fundamental properties used by each experiment.
    # Note: For Tokamaks and RFPs, property 'C' is considered 'True' for classification
    # purposes, as this is the only category available for current-driven devices (A).
    experiments = {
        "Tokamaks":               {"A": True, "B": True, "C": False},
        "LHD":                    {"A": False, "B": True, "C": True},
        "Wendelstein 7-X":        {"A": False, "B": True, "C": True},
        "NCSX":                   {"A": False, "B": True, "C": True},
        "Reversed Field Pinches": {"A": True, "B": True, "C": False}
    }

    # Lists to hold the names of experiments in each category
    only_b = []
    b_and_c = []
    a_b_and_c = []

    # Apply the classification logic
    for name, props in experiments.items():
        if props["A"]:
            # If property A (current) is present, it must fall into the "A, B, and C" category.
            a_b_and_c.append(name)
        elif props["B"] and props["C"]:
            # Devices using external 3D shaping without current are "B and C".
            b_and_c.append(name)
        elif props["B"] and not props["A"] and not props["C"]:
            # Devices using only shaping and a planar axis are "Only B".
            only_b.append(name)

    # Print the results in a clear format
    print("Classification of Fusion Experiments:")
    print("-" * 55)

    print("Category: Use only property B")
    if not only_b:
        print("  - None of the listed experiments")
    else:
        for exp in only_b:
            print(f"  - {exp}")
    print()

    print("Category: Use both properties B and C")
    if not b_and_c:
        print("  - None of the listed experiments")
    else:
        for exp in b_and_c:
            print(f"  - {exp}")
    print()

    print("Category: Use properties A, B, and C")
    if not a_b_and_c:
        print("  - None of the listed experiments")
    else:
        for exp in a_b_and_c:
            print(f"  - {exp}")
    print("-" * 55)

# Execute the function to print the classification
classify_fusion_experiments()