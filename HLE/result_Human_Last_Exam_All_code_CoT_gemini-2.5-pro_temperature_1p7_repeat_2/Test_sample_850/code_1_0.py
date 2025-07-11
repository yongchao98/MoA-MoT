import textwrap

def find_crystal_classes():
    """
    Analyzes the 32 crystal classes to find which ones are achiral, non-polar,
    and optically active, based on fundamental crystallographic principles.
    """
    # Data for the 32 point groups (crystal classes).
    # Format: (Hermann-Mauguin symbol, is_chiral, is_polar)
    # A crystal is optically active if and only if it is chiral.
    crystal_classes = [
        # Triclinic
        ("1", True, True), ("-1", False, False),
        # Monoclinic
        ("2", True, True), ("m", False, True), ("2/m", False, False),
        # Orthorhombic
        ("222", True, False), ("mm2", False, True), ("mmm", False, False),
        # Tetragonal
        ("4", True, True), ("-4", False, False), ("4/m", False, False),
        ("422", True, False), ("4mm", False, True), ("-42m", False, False),
        ("4/mmm", False, False),
        # Trigonal
        ("3", True, True), ("-3", False, False), ("32", True, False),
        ("3m", False, True), ("-3m", False, False),
        # Hexagonal
        ("6", True, True), ("-6", False, False), ("6/m", False, False),
        ("622", True, False), ("6mm", False, True), ("-6m2", False, False),
        ("6/mmm", False, False),
        # Cubic
        ("23", True, False), ("m-3", False, False), ("432", True, False),
        ("-43m", False, False), ("m-3m", False, False),
    ]

    # The condition for optical activity is chirality.
    # The user asks for classes that are:
    # 1. Achiral (is_chiral = False)
    # 2. Non-polar (is_polar = False)
    # 3. Optically Active (is_chiral = True)
    # This is a logical contradiction. The code will confirm this.

    achiral_nonpolar_classes = []
    optically_active_classes = []
    result_classes = []

    for name, is_chiral, is_polar in crystal_classes:
        is_optically_active = is_chiral

        if not is_chiral and not is_polar:
            achiral_nonpolar_classes.append(name)

        if is_optically_active:
            optically_active_classes.append(name)
        
        # Check for the contradictory condition
        if (not is_chiral) and (not is_polar) and is_optically_active:
            result_classes.append(name)

    print("--- Analysis of Crystal Classes for Optical Activity ---")
    print("\nPrinciple: A crystal class exhibits optical activity if and only if it is chiral.")
    print("Your query asks for classes that are simultaneously ACHIRAL and OPTICALLY ACTIVE, which is a contradiction.\n")

    print(f"First Group: Achiral and Non-polar crystal classes")
    print(f"   - Number of classes found: {len(achiral_nonpolar_classes)}")
    print("   - Members:", textwrap.fill(", ".join(achiral_nonpolar_classes), 70))
    
    print(f"\nSecond Group: Optically Active (i.e., Chiral) crystal classes")
    print(f"   - Number of classes found: {len(optically_active_classes)}")
    print("   - Members:", textwrap.fill(", ".join(optically_active_classes), 70))

    print("\n--- Final Result ---")
    print("The goal is to find the intersection of these two contradictory sets of requirements.")
    print("\nFinal Equation:")
    print(f"   Number of Achiral, Non-polar classes: {len(achiral_nonpolar_classes)}")
    print(f"   Number of Optically Active classes:  {len(optically_active_classes)}")
    print(f"   Number of classes meeting all criteria: {len(result_classes)}")
    print(f"\nConclusion: There are {len(result_classes)} achiral and non-polar crystal classes that are optically active.")

if __name__ == '__main__':
    find_crystal_classes()