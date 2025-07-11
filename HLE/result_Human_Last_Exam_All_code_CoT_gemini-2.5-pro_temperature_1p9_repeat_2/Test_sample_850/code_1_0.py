def find_contradictory_crystal_classes():
    """
    Explains the relationship between chirality, achirality, and optical activity in crystal classes.
    """
    # The 11 chiral crystal classes. These are the ONLY classes that can exhibit optical activity.
    # By definition, they lack mirror planes and inversion centers.
    chiral_optically_active_classes = ["1", "2", "3", "4", "6", "222", "32", "422", "622", "23", "432"]

    # The 21 achiral crystal classes. These classes are NOT optically active.
    # By definition, they all possess at least one mirror plane or an inversion center.
    achiral_non_optically_active_classes = [
        "-1", "2/m", "mmm",                               # Triclinic, Monoclinic, Orthorhombic
        "m", "mm2",                                       # (Monoclinic, Orthorhombic)
        "4/m", "4/mmm", "-4", "4mm", "-42m",             # Tetragonal
        "-3", "-3m",                                       # Trigonal
        "6/m", "6/mmm", "-6", "6mm", "-62m",             # Hexagonal
        "m-3", "m-3m"                                     # Cubic
    ]

    print("--- Understanding Crystal Symmetry and Optical Activity ---")
    print("\nStep 1: Define Optical Activity")
    print("For a crystal to be optically active, it must be CHIRAL.")
    print("A chiral crystal class has no center of inversion and no mirror planes.")
    print(f"The {len(chiral_optically_active_classes)} optically active (chiral) crystal classes are: {', '.join(chiral_optically_active_classes)}")

    print("\nStep 2: Define Achiral")
    print("An ACHIRAL crystal is the opposite of a chiral one.")
    print("An achiral crystal class MUST have a center of inversion or a mirror plane.")
    print(f"The {len(achiral_non_optically_active_classes)} achiral crystal classes are: {', '.join(achiral_non_optically_active_classes)}")

    print("\nStep 3: Analyze the User's Request")
    print("The request is for crystal classes that are both ACHIRAL and optically active.")

    print("\n--- Conclusion ---")
    print("A direct contradiction exists in the request.")
    print("A crystal CANNOT be both achiral and optically active simultaneously.")
    print("- Optical activity requires a CHIRAL structure.")
    print("- The request specifies an ACHIRAL structure.")
    print("\nTherefore, there are no crystal classes that are achiral but have the correct symmetry for optical activity.")
    print("The set of crystal classes that satisfy the request is empty.")

find_contradictory_crystal_classes()