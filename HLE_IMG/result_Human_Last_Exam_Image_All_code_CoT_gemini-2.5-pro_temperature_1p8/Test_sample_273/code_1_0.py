def determine_absolute_configuration():
    """
    This script explains the step-by-step process for determining the absolute
    configuration of the chiral centers in the provided molecule.
    """
    # Introduction
    print("To determine the absolute configuration of the molecule, we will apply the Cahn-Ingold-Prelog (CIP) rules to each of the three chiral centers.")
    print("We will label the chiral centers from left to right as C2, C3, and C4 for this explanation.\n")

    # --- Analysis for C2 ---
    print("--- Analysis of Chiral Center C2 (leftmost, with -OH) ---")
    print("1. Groups attached: -OH, -CH3, -H, and the rest of the molecule (R').")
    print("2. CIP Priority Assignment:")
    print("   - Priority 1: -OH (Oxygen has the highest atomic number)")
    print("   - Priority 2: R' (The rest of the molecule, -C(Et)-...)")
    print("   - Priority 3: -CH3")
    print("   - Priority 4: -H (Hydrogen has the lowest atomic number)")
    print("3. 3D Orientation:")
    print("   - In the drawing, the -OH and -CH3 groups have dashed bonds (pointing away). For a tetrahedral carbon, this implies the unseen -H atom must have a wedged bond (pointing towards the viewer).")
    print("   - Therefore, the lowest priority group (-H) is pointing towards the viewer.")
    print("4. Determining Configuration:")
    print("   - The path from priority 1 -> 2 -> 3 (-OH -> R' -> -CH3) is counter-clockwise.")
    print("   - Rule: When the lowest priority group points towards the viewer, we reverse the assignment. Counter-clockwise (S) becomes R.")
    print("=> Configuration of C2 is (R).\n")

    # --- Analysis for C3 ---
    print("--- Analysis of Chiral Center C3 (middle, with ethyl group) ---")
    print("1. Groups attached: -CH2CH3 (ethyl), -H, the group to the left (-CH(OH)CH3), and the group to the right (-CH(Me)CH2NH2).")
    print("2. CIP Priority Assignment:")
    print("   - Priority 1: Left group (-CH(OH)CH3, carbon bonded to O, C, H)")
    print("   - Priority 2: Right group (-CH(Me)CH2NH2, carbon bonded to C, C, H)")
    print("   - Priority 3: Ethyl group (-CH2CH3, carbon bonded to C, H, H)")
    print("   - Priority 4: -H")
    print("3. 3D Orientation:")
    print("   - The ethyl group has a wedged bond (towards the viewer). This implies the unseen -H atom has a dashed bond (away from the viewer).")
    print("   - The lowest priority group (-H) is pointing away from the viewer.")
    print("4. Determining Configuration:")
    print("   - The path from priority 1 -> 2 -> 3 (Left group -> Right group -> Ethyl group) is clockwise.")
    print("   - Rule: When the lowest priority group points away, the assignment is direct. Clockwise is R.")
    print("=> Configuration of C3 is (R).\n")

    # --- Analysis for C4 ---
    print("--- Analysis of Chiral Center C4 (rightmost) ---")
    print("1. Groups attached: -CH2NH2, -CH3, -H, and the rest of the molecule (R'').")
    print("2. CIP Priority Assignment:")
    print("   - Priority 1: -CH2NH2 (carbon bonded to N)")
    print("   - Priority 2: R'' (the rest of the molecule, -CH(Et)-...)")
    print("   - Priority 3: -CH3")
    print("   - Priority 4: -H")
    print("3. 3D Orientation:")
    print("   - The -CH3 group has a dashed bond (away from viewer). This implies the unseen -H atom has a wedged bond (towards the viewer).")
    print("   - The lowest priority group (-H) is pointing towards the viewer.")
    print("4. Determining Configuration:")
    print("   - The path from priority 1 -> 2 -> 3 (-CH2NH2 -> R'' -> -CH3) is counter-clockwise.")
    print("   - Rule: We reverse the assignment. Counter-clockwise (S) becomes R.")
    print("=> Configuration of C4 is (R).\n")

    # --- Final Result ---
    print("--- Final Configuration ---")
    print("The final configuration is the combination of the individual configurations at centers 2, 3, and 4:")
    print("Configuration at Center 2 = R")
    print("Configuration at Center 3 = R")
    print("Configuration at Center 4 = R")

# Run the analysis
determine_absolute_configuration()