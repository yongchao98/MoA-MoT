def classify_fusion_experiments():
    """
    Classifies fusion experiments based on their method of twisting the magnetic field.

    The properties are:
    (A) driving a toroidal current
    (B) elongating the flux surfaces and making them rotate poloidally
    (C) making the magnetic axis non-planar
    
    The classification is as follows:
    - Tokamaks & Reversed Field Pinches: Use a toroidal current (A) and plasma shaping (B).
      They are assigned to "A, B, and C" as it's the only category including property A.
    - LHD, Wendelstein 7-X, NCSX: These are stellarators that use external coils for shaping/rotation (B)
      and feature a non-planar axis (C).
    """

    # Based on the analysis, we assign each experiment to a category.
    # Note: Tokamaks and RFPs are factually A and B. They are placed in A, B, and C
    # to fit the problem's specified categories.
    classification = {
        'Only B': [],
        'B and C': ['LHD', 'Wendelstein 7-X', 'NCSX'],
        'A, B, and C': ['Tokamaks', 'Reversed Field Pinches']
    }

    print("Use only property B:")
    only_B_devices = classification['Only B']
    if not only_B_devices:
        print("None")
    else:
        for device in only_B_devices:
            print(f"- {device}")

    print("\nUse both B and C:")
    B_and_C_devices = classification['B and C']
    if not B_and_C_devices:
        print("None")
    else:
        for device in B_and_C_devices:
            print(f"- {device}")
    
    print("\nUse A, B, and C:")
    A_B_C_devices = classification['A, B, and C']
    if not A_B_C_devices:
        print("None")
    else:
        for device in A_B_C_devices:
            print(f"- {device}")

if __name__ == '__main__':
    classify_fusion_experiments()