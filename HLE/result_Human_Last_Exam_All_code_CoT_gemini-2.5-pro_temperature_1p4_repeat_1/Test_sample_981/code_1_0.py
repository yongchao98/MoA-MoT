def classify_fusion_experiments():
    """
    Classifies fusion experiments based on how they twist their magnetic fields.

    The properties are:
    (A) Driving a toroidal current.
    (B) Elongating flux surfaces and making them rotate poloidally.
    (C) Making the magnetic axis non-planar.

    This code implements the logic that:
    - Stellarators (LHD, W7-X, NCSX) use 3D geometry (B and C) but not a large driven current.
    - Tokamaks use current (A), geometric shaping (B), and 3D fields for control (C).
    - RFPs use current (A) and rely on an essential 3D dynamo effect (interpreted as B and C).
    """

    # Represent each experiment and its properties
    # This dictionary encodes the analysis explained in the plan.
    devices = {
        "Tokamaks":               {"A": True, "B": True, "C": True},
        "LHD":                    {"A": False, "B": True, "C": True},
        "Wendelstein 7-X":        {"A": False, "B": True, "C": True},
        "NCSX":                   {"A": False, "B": True, "C": True},
        "Reversed Field Pinches": {"A": True, "B": True, "C": True}
    }

    # Initialize lists for each category
    only_B = []
    B_and_C = []
    A_B_and_C = []

    # Categorize each device based on its properties
    for name, props in devices.items():
        use_A = props["A"]
        use_B = props["B"]
        use_C = props["C"]

        # Check for category "A, B, and C"
        if use_A and use_B and use_C:
            A_B_and_C.append(name)
        # Check for category "B and C"
        elif not use_A and use_B and use_C:
            B_and_C.append(name)
        # Check for category "Only B"
        elif not use_A and use_B and not use_C:
            only_B.append(name)

    # Print the final classification
    print("Classification of Fusion Experiments by Magnetic Field Properties:")
    print("="*60)

    print("\nExperiments that use only property B (Elongated/rotating flux surfaces):")
    if not only_B:
        print("  - None")
    else:
        for device in only_B:
            print(f"  - {device}")

    print("\nExperiments that use both properties B and C (Non-planar axis):")
    if not B_and_C:
        print("  - None")
    else:
        for device in B_and_C:
            print(f"  - {device}")

    print("\nExperiments that use properties A (Current), B, and C:")
    if not A_B_and_C:
        print("  - None")
    else:
        for device in A_B_and_C:
            print(f"  - {device}")

if __name__ == '__main__':
    classify_fusion_experiments()