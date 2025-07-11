def classify_fusion_experiments():
    """
    Classifies fusion experiments based on their method of twisting the magnetic field.

    The properties are:
    (A) driving a toroidal current
    (B) elongating the flux surfaces and making them rotate poloidally
    (C) making the magnetic axis non-planar
    
    The classification is based on the following reasoning:
    - 'Only property B': Devices with a planar axis where the defining characteristic is the
      rotation of flux surfaces, regardless of the method (e.g., current-driven). This
      group includes Tokamaks and RFPs.
    - 'Both B and C': Advanced stellarators that use 3D shaping and a non-planar axis
      to generate the twist, without relying on a major driven current. This group
      includes LHD and Wendelstein 7-X.
    - 'A, B, and C': Hybrid or advanced concepts designed to operate with significant
      plasma current in a non-planar, 3D-shaped magnetic field. NCSX fits this as
      a compact stellarator concept designed for high bootstrap current.
    """
    
    classifications = {
        "Only property B": ["Tokamaks", "Reversed Field Pinches"],
        "Both B and C": ["LHD", "Wendelstein 7-X"],
        "A, B, and C": ["NCSX"]
    }
    
    print("Classification of Fusion Experiments:")
    print("-" * 35)

    for category, experiments in classifications.items():
        print(f"\nCategory: '{category}'")
        for exp in experiments:
            print(f" - {exp}")

classify_fusion_experiments()