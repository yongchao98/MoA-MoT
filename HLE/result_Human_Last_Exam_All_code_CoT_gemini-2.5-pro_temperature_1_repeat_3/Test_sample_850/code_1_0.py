def analyze_crystal_classes():
    """
    Analyzes the 32 crystallographic point groups to find which are achiral, 
    non-polar, and optically active, demonstrating the inherent contradiction in the request.
    """
    # Data for the 32 crystallographic point groups.
    # 'is_chiral': True if the group lacks mirror planes and inversion centers.
    #              This is the sole condition for optical activity.
    # 'is_polar': True if the group has a unique polar axis (pyroelectric).
    point_groups = [
        # Triclinic
        {'name': '1', 'is_chiral': True, 'is_polar': True},
        {'name': '-1', 'is_chiral': False, 'is_polar': False},
        # Monoclinic
        {'name': '2', 'is_chiral': True, 'is_polar': True},
        {'name': 'm', 'is_chiral': False, 'is_polar': True},
        {'name': '2/m', 'is_chiral': False, 'is_polar': False},
        # Orthorhombic
        {'name': '222', 'is_chiral': True, 'is_polar': False},
        {'name': 'mm2', 'is_chiral': False, 'is_polar': True},
        {'name': 'mmm', 'is_chiral': False, 'is_polar': False},
        # Tetragonal
        {'name': '4', 'is_chiral': True, 'is_polar': True},
        {'name': '-4', 'is_chiral': False, 'is_polar': False},
        {'name': '4/m', 'is_chiral': False, 'is_polar': False},
        {'name': '422', 'is_chiral': True, 'is_polar': False},
        {'name': '4mm', 'is_chiral': False, 'is_polar': True},
        {'name': '-42m', 'is_chiral': False, 'is_polar': False},
        {'name': '4/mmm', 'is_chiral': False, 'is_polar': False},
        # Trigonal
        {'name': '3', 'is_chiral': True, 'is_polar': True},
        {'name': '-3', 'is_chiral': False, 'is_polar': False},
        {'name': '32', 'is_chiral': True, 'is_polar': False},
        {'name': '3m', 'is_chiral': False, 'is_polar': True},
        {'name': '-3m', 'is_chiral': False, 'is_polar': False},
        # Hexagonal
        {'name': '6', 'is_chiral': True, 'is_polar': True},
        {'name': '-6', 'is_chiral': False, 'is_polar': False},
        {'name': '6/m', 'is_chiral': False, 'is_polar': False},
        {'name': '622', 'is_chiral': True, 'is_polar': False},
        {'name': '6mm', 'is_chiral': False, 'is_polar': True},
        {'name': '-6m2', 'is_chiral': False, 'is_polar': False}, # same as -62m
        {'name': '6/mmm', 'is_chiral': False, 'is_polar': False},
        # Cubic
        {'name': '23', 'is_chiral': True, 'is_polar': False},
        {'name': 'm-3', 'is_chiral': False, 'is_polar': False},
        {'name': '432', 'is_chiral': True, 'is_polar': False},
        {'name': '-43m', 'is_chiral': False, 'is_polar': False},
        {'name': 'm-3m', 'is_chiral': False, 'is_polar': False},
    ]

    # --- Step 1: Explain the physics ---
    print("--- Analysis of Crystal Symmetry and Optical Activity ---")
    print("A crystal class can exhibit optical activity if and only if it is chiral.")
    print("A chiral class is one that lacks a center of inversion and mirror planes.")
    print("The user's question asks for classes that are both 'achiral' and 'optically active', which is a contradiction.\n")

    # --- Step 2: Find classes meeting the user's criteria (achiral and non-polar) ---
    achiral_non_polar_classes = [
        pg['name'] for pg in point_groups if not pg['is_chiral'] and not pg['is_polar']
    ]
    print(f"Step A: Finding crystal classes that are ACHIRAL and NON-POLAR...")
    print(f"Result ({len(achiral_non_polar_classes)} classes found): {', '.join(achiral_non_polar_classes)}\n")

    # --- Step 3: Find classes that are optically active (chiral) ---
    optically_active_classes = [pg['name'] for pg in point_groups if pg['is_chiral']]
    print(f"Step B: Identifying all OPTICALLY ACTIVE (i.e., CHIRAL) crystal classes...")
    print(f"Result ({len(optically_active_classes)} classes found): {', '.join(optically_active_classes)}\n")

    # --- Step 4: Show the intersection is empty ---
    print("Step C: Checking for classes that meet both conditions (from Step A and Step B)...")
    intersection = set(achiral_non_polar_classes).intersection(set(optically_active_classes))
    
    if not intersection:
        print("Result: The intersection is empty.")
        print("\n--- Conclusion ---")
        print("There are no crystal classes that are both achiral and optically active.")
        print("The conditions are mutually exclusive.")
    else:
        # This case is physically impossible and will not be reached.
        print(f"Found impossible classes: {', '.join(intersection)}")

if __name__ == '__main__':
    analyze_crystal_classes()