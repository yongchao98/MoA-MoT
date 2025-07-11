def determine_absolute_configuration():
    """
    This script provides a step-by-step determination of the absolute
    configuration (R/S) for each chiral center in the provided molecule.
    """
    print("Determining the absolute configuration of the molecule's three chiral centers.")
    print("The centers are analyzed from left to right.\n")

    # --- Center 1: Carbon with -OH ---
    print("--- Analysis of Center 1 (leftmost, with -OH) ---")
    print("Groups: -OH, -CH3, -H (implicit), and the rest of the molecule (R-group).")
    print("CIP Priorities:")
    print("1: -OH (Oxygen has the highest atomic number)")
    print("2: R-group (the rest of the molecule, connected via Carbon)")
    print("3: -CH3 (methyl group)")
    print("4: -H (Hydrogen has the lowest atomic number)")
    print("\nGeometry and Configuration:")
    print("The -OH (1) is on a wedge, and -CH3 (3) is on a dash. The -H (4) is in the plane.")
    print("To determine the configuration, we can swap the lowest priority group (-H, #4) with the group on the dash (-CH3, #3).")
    print("After the swap, -H is in the back. The sequence 1 -> 2 -> 3 is now counter-clockwise (S).")
    print("Since we performed one swap, we reverse the result.")
    print("Result for Center 1: S becomes R. The configuration is R.\n")

    # --- Center 2: Carbon with Ethyl ---
    print("--- Analysis of Center 2 (middle, with ethyl group) ---")
    print("Groups: Ethyl (-CH2CH3), left R-group (-CH(OH)CH3), right R-group (-CH(CH3)CH2NH2), and -H (implicit).")
    print("CIP Priorities:")
    print("1: Left R-group (-CH(OH)CH3) (due to the Oxygen atom)")
    print("2: Right R-group (-CH(CH3)CH2NH2) (due to Nitrogen vs Carbon in ethyl)")
    print("3: Ethyl group (-CH2CH3)")
    print("4: -H")
    print("\nGeometry and Configuration:")
    print("The ethyl group (3) is on a wedge, and the lowest priority group -H (4) is on a dash (in the back).")
    print("Since the lowest priority group is in the back, we trace the path 1 -> 2 -> 3.")
    print("The path is clockwise.")
    print("Result for Center 2: Clockwise with #4 in the back is R. The configuration is R.\n")

    # --- Center 3: Carbon with -CH2NH2 ---
    print("--- Analysis of Center 3 (rightmost) ---")
    print("Groups: -CH2NH2, -CH3, left R-group, and -H (implicit).")
    print("CIP Priorities:")
    print("1: -CH2NH2 (due to the Nitrogen atom)")
    print("2: Left R-group")
    print("3: -CH3")
    print("4: -H")
    print("\nGeometry and Configuration:")
    print("The -CH3 (3) is on a dash, which means the implicit -H (4) must be on a wedge (pointing out).")
    print("We trace the path 1 -> 2 -> 3, which is counter-clockwise.")
    print("Since the lowest priority group (-H) is on a wedge, we must reverse the result.")
    print("Result for Center 3: Counter-clockwise would be S, so reversing it gives R. The configuration is R.\n")

    # --- Final Conclusion ---
    print("--- Final Conclusion ---")
    print("The absolute configurations for the three chiral centers, from left to right, are:")
    print("Center 1: R")
    print("Center 2: R")
    print("Center 3: R")

determine_absolute_configuration()
<<<R,R,R>>>