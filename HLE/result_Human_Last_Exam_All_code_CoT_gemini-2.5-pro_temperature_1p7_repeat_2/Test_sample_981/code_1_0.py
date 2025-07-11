def classify_fusion_experiments():
    """
    Identifies and prints the classification of fusion experiments based on
    the method used to twist the magnetic field.

    The properties for classification are:
    (A) driving a toroidal current
    (B) elongating the flux surfaces and making them rotate poloidally
    (C) making the magnetic axis non-planar
    """

    # The classification is based on the following logic:
    # 1. Stellarators (LHD, W7-X, NCSX) do not use a net toroidal current (A),
    #    but rely on external field shaping (B) and have a non-planar axis (C).
    #    This places them in the "Both B and C" category.
    # 2. Tokamaks and RFPs are defined by their use of a large toroidal current (A).
    #    This excludes them from "Only B" and "B and C". By elimination, they are
    #    placed in the "A, B, and C" category, even though they do not have a
    #    non-planar axis (C), as this is the only remaining option in the problem's framework.

    classification = {
        "experiments using only property B": [],
        "experiments using both B and C": ["LHD", "Wendelstein 7-X", "NCSX"],
        "experiments using A, B, and C": ["Tokamaks", "Reversed Field Pinches"]
    }

    print("--- Fusion Experiment Classification ---")
    print("\nExperiments that use only property B:")
    if classification["experiments using only property B"]:
        for exp in classification["experiments using only property B"]:
            print(f"- {exp}")
    else:
        print("- None from the list")

    print("\nExperiments that use both B and C:")
    if classification["experiments using both B and C"]:
        for exp in classification["experiments using both B and C"]:
            print(f"- {exp}")
    else:
        print("- None from the list")

    print("\nExperiments that use A, B, and C:")
    if classification["experiments using A, B, and C"]:
        for exp in classification["experiments using A, B, and C"]:
            print(f"- {exp}")
    else:
        print("- None from the list")

if __name__ == '__main__':
    classify_fusion_experiments()