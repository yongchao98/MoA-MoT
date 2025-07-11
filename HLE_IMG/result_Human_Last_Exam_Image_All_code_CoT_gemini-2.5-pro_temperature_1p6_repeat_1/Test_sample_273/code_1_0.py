def determine_configuration():
    """
    This script determines the absolute configuration (R/S) for the chiral centers
    of the given molecule.
    """
    print("Step 1: IUPAC Name and Chiral Center Identification")
    # Based on IUPAC rules, the principal functional group is the alcohol (-OH).
    # The longest carbon chain containing the C-OH carbon is a 5-carbon pentane chain.
    # Numbering from left to right gives the -OH group the lowest number (2).
    # The substituents are an ethyl group at C3 and an aminomethyl group at C4.
    # The image corresponds to a stereoisomer of 4-(aminomethyl)-3-ethylpentan-2-ol.
    # Let's refer to the carbons by their number in this chain.
    iupac_name = "4-(aminomethyl)-3-ethylpentan-2-ol"
    chiral_centers = ["C2", "C3", "C4"]
    print(f"The molecule is identified as an isomer of: {iupac_name}")
    print(f"The chiral centers are at positions: {', '.join(chiral_centers)}\n")

    # --- Analysis for C2 ---
    print("--- Analysis for Chiral Center C2 ---")
    print("Substituents on C2: -OH, -CH3 (C1), -C3_group, -H")
    print("Priority Assignment (CIP Rules):")
    print("1: -OH (highest atomic number O vs C)")
    print("2: -C3_group (-CH(Et)-...)")
    print("3: -CH3 (C1)")
    print("4: -H (lowest atomic number)")
    print("Spatial Orientation: The lowest priority group (-H) is pointing forwards (wedged).")
    print("Path from priority 1 -> 2 -> 3 is counter-clockwise.")
    print("Result: Counter-clockwise with H in front gives R configuration.")
    c2_config = "2R"
    print(f"Configuration at C2 is: {c2_config}\n")

    # --- Analysis for C3 ---
    print("--- Analysis for Chiral Center C3 ---")
    print("Substituents on C3: -C2_group, -C4_group, -Ethyl, -H")
    print("Priority Assignment (CIP Rules):")
    print("1: -C2_group (-CH(OH)-... contains O)")
    print("2: -C4_group (-CH(CH2NH2)-...)")
    print("3: -Ethyl (-CH2CH3)")
    print("4: -H")
    print("Spatial Orientation: The lowest priority group (-H) is pointing backwards (dashed).")
    print("Path from priority 1 -> 2 -> 3 is clockwise.")
    print("Result: Clockwise with H in back gives R configuration.")
    c3_config = "3R"
    print(f"Configuration at C3 is: {c3_config}\n")

    # --- Analysis for C4 ---
    print("--- Analysis for Chiral Center C4 ---")
    print("Substituents on C4: -CH2NH2, -C3_group, -CH3, -H")
    print("Priority Assignment (CIP Rules):")
    print("1: -CH2NH2 (C is bonded to N)")
    print("2: -C3_group (C is bonded to C,C,H)")
    print("3: -CH3 (C is bonded to H,H,H)")
    print("4: -H")
    print("Spatial Orientation: The lowest priority group (-H) is pointing forwards (wedged).")
    print("Path from priority 1 -> 2 -> 3 is clockwise.")
    print("Result: Clockwise with H in front gives S configuration.")
    c4_config = "4S"
    print(f"Configuration at C4 is: {c4_config}\n")
    
    # --- Final Conclusion ---
    print("--- Final Configuration ---")
    final_config = f"({c2_config}, {c3_config}, {c4_config})"
    print(f"The absolute configuration of the molecule is {final_config}.")

if __name__ == '__main__':
    determine_configuration()