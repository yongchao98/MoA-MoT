def determine_absolute_configuration():
    """
    This script determines and explains the absolute configuration of the given molecule.
    """
    print("Step 1: IUPAC Naming and Identification of Chiral Centers")
    print("The principal functional group is the alcohol (-OH). The longest carbon chain containing the C-OH group is a pentane chain.")
    print("Numbering from left to right gives the -OH group the lowest number (2).")
    print("The IUPAC name is: 4-(aminomethyl)-3-ethylpentan-2-ol.")
    print("The molecule has three chiral centers at carbons C2, C3, and C4.\n")

    print("Step 2: Analysis of Each Chiral Center (Assuming -OH on C2 is Wedged)")
    print("Note: The original drawing for C2 is ambiguous. We assume the -OH is wedged and -CH3 is dashed.\n")

    # --- Analysis of C2 ---
    print("--- Configuration at C2 ---")
    print("Attached groups: -OH, -C3, -CH3(C1), -H")
    print("Priorities (CIP Rules):")
    print("1: -OH")
    print("2: -C3")
    print("3: -CH3")
    print("4: -H")
    print("Arrangement: -OH(1) is wedged, -CH3(3) is dashed. The implied -H(4) is in the plane.")
    print("Viewing the molecule with the lowest priority group (-H) pointing away, the path from 1 -> 2 -> 3 is counter-clockwise.")
    print("Therefore, the configuration at C2 is S.\n")

    # --- Analysis of C3 ---
    print("--- Configuration at C3 ---")
    print("Attached groups: -C2, -C4, -Ethyl, -H")
    print("Priorities (CIP Rules):")
    print("1: -C2 (attached to O)")
    print("2: -C4 (attached to C, C, H where one path leads to N)")
    print("3: -Ethyl (attached to C, H, H)")
    print("4: -H")
    print("Arrangement: Ethyl(3) is wedged, so the implied -H(4) is dashed (away).")
    print("With the lowest priority group away, the path from 1 -> 2 -> 3 is counter-clockwise.")
    print("Therefore, the configuration at C3 is S.\n")

    # --- Analysis of C4 ---
    print("--- Configuration at C4 ---")
    print("Attached groups: -CH2NH2, -C3, -CH3(C5), -H")
    print("Priorities (CIP Rules):")
    print("1: -CH2NH2 (C attached to N)")
    print("2: -C3")
    print("3: -CH3")
    print("4: -H")
    print("Arrangement: -CH3(3) is dashed, so the implied -H(4) is wedged (forward).")
    print("The path from 1 -> 2 -> 3 is clockwise.")
    print("Since the lowest priority group (-H) is forward, we reverse the result (clockwise=R becomes S).")
    print("Therefore, the configuration at C4 is S.\n")

    print("--- Final Result ---")
    c2_num, c2_config = 2, 'S'
    c3_num, c3_config = 3, 'S'
    c4_num, c4_config = 4, 'S'
    
    print(f"The final absolute configuration of the molecule is ({c2_num}{c2_config}, {c3_num}{c3_config}, {c4_num}{c4_config}).")

determine_absolute_configuration()
<<<(2S, 3S, 4S)>>>