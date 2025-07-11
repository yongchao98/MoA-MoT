def solve_stereochemistry():
    """
    This script determines the absolute configuration of the provided molecule
    by analyzing each chiral center according to Cahn-Ingold-Prelog (CIP) rules.
    """
    print("Determining the absolute configuration of the molecule...")
    print("The molecule has three chiral centers, which we label C2, C3, and C4 from left to right.")
    print("-" * 40)

    # --- Analysis of Chiral Center C2 ---
    print("Analysis of Chiral Center C2 (with -OH group):")
    print("Substituent Priorities (CIP Rules):")
    print("1: -OH")
    print("2: -C3H(Et)... (rest of chain)")
    print("3: -CH3")
    print("4: -H (implicit)")
    print("\n3D Arrangement:")
    print("-OH is dashed (away), so -H is wedged (forward).")
    print("Path from 1 -> 2 -> 3 is clockwise (R).")
    print("Since the lowest priority group (4) is wedged, we reverse R to S.")
    c2_config = "S"
    print(f"Result for C2: ({c2_config})")
    print("-" * 40)

    # --- Analysis of Chiral Center C3 ---
    print("Analysis of Chiral Center C3 (with ethyl group):")
    print("Substituent Priorities (CIP Rules):")
    print("1: -C2H(OH)... (group with hydroxyl)")
    print("2: -C4H(Me)... (group with amine)")
    print("3: -CH2CH3 (Ethyl)")
    print("4: -H (implicit)")
    print("\n3D Arrangement:")
    print("Ethyl group is wedged (forward), so -H is dashed (away).")
    print("Path from 1 -> 2 -> 3 is counter-clockwise (S).")
    print("Since the lowest priority group (4) is dashed, the configuration is S.")
    c3_config = "S"
    print(f"Result for C3: ({c3_config})")
    print("-" * 40)

    # --- Analysis of Chiral Center C4 ---
    print("Analysis of Chiral Center C4 (with methyl and aminoethyl groups):")
    print("Substituent Priorities (CIP Rules):")
    print("1: -CH2NH2 (C bonded to N gives highest priority)")
    print("2: -C3H(Et)... (rest of chain)")
    print("3: -CH3")
    print("4: -H (implicit)")
    print("\n3D Arrangement:")
    print("-CH3 is dashed (away), so -H is wedged (forward).")
    print("Path from 1 -> 2 -> 3 is counter-clockwise (S).")
    print("Since the lowest priority group (4) is wedged, we reverse S to R.")
    c4_config = "R"
    print(f"Result for C4: ({c4_config})")
    print("-" * 40)

    # --- Final Answer ---
    final_config_string = f"(2{c2_config}, 3{c3_config}, 4{c4_config})"
    print(f"The complete absolute configuration of the molecule is: {final_config_string}")

solve_stereochemistry()