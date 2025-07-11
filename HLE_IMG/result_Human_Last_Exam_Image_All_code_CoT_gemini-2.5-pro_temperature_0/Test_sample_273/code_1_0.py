def determine_absolute_configuration():
    """
    This script outlines the step-by-step process to determine the absolute
    configuration of the chiral centers in the provided molecule.
    """
    print("Determining the absolute configuration of the molecule:")
    print("The molecule has three chiral centers. We will analyze them from left to right (C2, C3, C4).\n")

    # --- Analysis of C2 (leftmost chiral center) ---
    print("--- Chiral Center C2 (Carbon with -OH and -CH3) ---")
    print("1. Priorities:")
    print("   1: -OH")
    print("   2: -C3 (rest of the chain to the right)")
    print("   3: -CH3")
    print("   4: -H (implicit)")
    print("2. 3D Orientation:")
    print("   - The -OH group is on a dash (away).")
    print("   - The implicit -H atom is on a wedge (forward).")
    print("3. Assignment:")
    print("   - With the lowest priority group (-H) pointing forward, the path 1 -> 2 -> 3 is counter-clockwise (S).")
    print("   - Since the H is forward, we reverse the result.")
    print("   - Configuration at C2 is R.\n")

    # --- Analysis of C3 (middle chiral center) ---
    print("--- Chiral Center C3 (Carbon with ethyl group) ---")
    print("1. Priorities:")
    print("   1: -C2 (carbon attached to O)")
    print("   2: -C4 (carbon attached to N via another C)")
    print("   3: -CH2CH3 (ethyl)")
    print("   4: -H (implicit)")
    print("2. 3D Orientation:")
    print("   - The ethyl group is on a wedge (forward).")
    print("   - The implicit -H atom is on a dash (away).")
    print("3. Assignment:")
    print("   - With the lowest priority group (-H) pointing away, the path 1 -> 2 -> 3 is clockwise (R).")
    print("   - Since the H is away, the result is as observed.")
    print("   - Configuration at C3 is R.\n")

    # --- Analysis of C4 (rightmost chiral center) ---
    print("--- Chiral Center C4 (Carbon with -CH3 and -CH2NH2) ---")
    print("1. Priorities:")
    print("   1: -CH2NH2 (carbon attached to N)")
    print("   2: -C3 (rest of the chain to the left)")
    print("   3: -CH3")
    print("   4: -H (implicit)")
    print("2. 3D Orientation:")
    print("   - The -CH3 group is on a dash (away).")
    print("   - The implicit -H atom is on a wedge (forward).")
    print("3. Assignment:")
    print("   - With the lowest priority group (-H) pointing forward, the path 1 -> 2 -> 3 is counter-clockwise (S).")
    print("   - Since the H is forward, we reverse the result.")
    print("   - Configuration at C4 is R.\n")

    print("="*30)
    print("Final Conclusion:")
    print("The absolute configuration of the molecule is (2R, 3R, 4R).")
    print("="*30)

determine_absolute_configuration()
<<<2R, 3R, 4R>>>