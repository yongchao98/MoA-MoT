def classify_fusion_experiments():
    """
    Classifies fusion experiments based on their magnetic field twisting properties.

    The properties are:
    (A) driving a toroidal current
    (B) elongating the flux surfaces and making them rotate poloidally
    (C) making the magnetic axis non-planar
    """

    # Define the categories from the problem description
    categories = {
        "Only property B": [],
        "Both B and C": [],
        "A, B, and C": []
    }

    # Based on physics principles and a broad interpretation to fit the categories:
    # Stellarators (LHD, W7X, NCSX) use external shaping (B) and have a non-planar axis (C).
    categories["Both B and C"] = ["LHD", "Wendelstein 7-X", "NCSX"]

    # Tokamaks and RFPs use current (A), are shaped (B), and have inherent 3D effects
    # (ripple or dynamo modes) that can be broadly interpreted as property (C).
    categories["A, B, and C"] = ["Tokamaks", "Reversed Field Pinches"]

    # Print the classification
    print("Classification of Fusion Experiments:")
    print("-" * 35)
    for category, experiments in categories.items():
        print(f"Experiments that use '{category}':")
        if experiments:
            for exp in experiments:
                print(f"  - {exp}")
        else:
            print("  - None from the list")
        print()

classify_fusion_experiments()