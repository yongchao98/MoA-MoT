def find_crystal_classes():
    """
    Analyzes the 32 crystallographic point groups to find those that meet
    the user's criteria and explains the result.
    """
    # Data for the 32 point groups (crystal classes)
    # is_chiral: True if the point group lacks mirror planes and inversion centers.
    #            This is the necessary condition for optical activity.
    # is_polar: True if the point group has a net dipole moment.
    point_groups = [
        # Triclinic
        {"name": "1", "is_chiral": True, "is_polar": True},
        {"name": "-1", "is_chiral": False, "is_polar": False},
        # Monoclinic
        {"name": "2", "is_chiral": True, "is_polar": True},
        {"name": "m", "is_chiral": False, "is_polar": True},
        {"name": "2/m", "is_chiral": False, "is_polar": False},
        # Orthorhombic
        {"name": "222", "is_chiral": True, "is_polar": False},
        {"name": "mm2", "is_chiral": False, "is_polar": True},
        {"name": "mmm", "is_chiral": False, "is_polar": False},
        # Tetragonal
        {"name": "4", "is_chiral": True, "is_polar": True},
        {"name": "-4", "is_chiral": False, "is_polar": False},
        {"name": "4/m", "is_chiral": False, "is_polar": False},
        {"name": "422", "is_chiral": True, "is_polar": False},
        {"name": "4mm", "is_chiral": False, "is_polar": True},
        {"name": "-42m", "is_chiral": False, "is_polar": False},
        {"name": "4/mmm", "is_chiral": False, "is_polar": False},
        # Trigonal
        {"name": "3", "is_chiral": True, "is_polar": True},
        {"name": "-3", "is_chiral": False, "is_polar": False},
        {"name": "32", "is_chiral": True, "is_polar": False},
        {"name": "3m", "is_chiral": False, "is_polar": True},
        {"name": "-3m", "is_chiral": False, "is_polar": False},
        # Hexagonal
        {"name": "6", "is_chiral": True, "is_polar": True},
        {"name": "-6", "is_chiral": False, "is_polar": False},
        {"name": "6/m", "is_chiral": False, "is_polar": False},
        {"name": "622", "is_chiral": True, "is_polar": False},
        {"name": "6mm", "is_chiral": False, "is_polar": True},
        {"name": "-6m2", "is_chiral": False, "is_polar": False},
        {"name": "6/mmm", "is_chiral": False, "is_polar": False},
        # Cubic
        {"name": "23", "is_chiral": True, "is_polar": False},
        {"name": "m-3", "is_chiral": False, "is_polar": False},
        {"name": "432", "is_chiral": True, "is_polar": False},
        {"name": "-43m", "is_chiral": False, "is_polar": False},
        {"name": "m-3m", "is_chiral": False, "is_polar": False},
    ]

    # The condition for optical activity is chirality.
    # The user asks for classes that are:
    # 1. Achiral (is_chiral = False)
    # 2. Non-polar (is_polar = False)
    # 3. Optically active (is_chiral = True)
    # The combination of (1) and (3) is a contradiction.

    print("Analyzing the request: Find crystal classes that are achiral, non-polar, and optically active.\n")
    print("Step 1: The condition for a crystal class to be optically active is that it must be chiral.")
    print("Step 2: The condition for a class to be chiral is that it must lack all improper rotation symmetries (no mirror planes or centers of inversion).")
    print("Step 3: The user's request asks for classes that are simultaneously achiral (possess improper symmetry) and optically active (lack improper symmetry).")
    print("\nThis is a fundamental contradiction. A crystal class cannot be both achiral and optically active.\n")

    # Let's find the classes that meet the first two criteria (achiral and non-polar)
    # to show why they cannot be optically active.
    achiral_nonpolar_classes = [
        pg for pg in point_groups if not pg["is_chiral"] and not pg["is_polar"]
    ]

    print("The following crystal classes are achiral and non-polar.")
    print("Because they are achiral, they CANNOT be optically active:")
    class_names = [pg['name'] for pg in achiral_nonpolar_classes]
    print(', '.join(class_names))

    print("\nConclusion: There are 0 crystal classes that satisfy all the requested conditions.")

if __name__ == '__main__':
    find_crystal_classes()