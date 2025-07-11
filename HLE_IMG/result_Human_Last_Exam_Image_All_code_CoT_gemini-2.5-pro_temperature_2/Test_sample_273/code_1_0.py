def solve_stereochemistry():
    """
    This script determines the absolute configuration (R/S) of the three stereocenters
    in the provided molecule by applying the Cahn-Ingold-Prelog (CIP) priority rules.
    The stereocenters are labeled C2, C3, and C4 from left to right for clarity.
    """

    print("Analyzing the stereocenter at C2 (leftmost, with -OH group):")
    print("Step 1: Assign CIP priorities based on atomic number at the first point of difference.")
    print("1: -OH (Oxygen Z=8)")
    print("2: -C(Et)... (the rest of the chain, a carbon bonded to two other carbons)")
    print("3: -CH3 (a carbon bonded to three hydrogens)")
    print("4: -H (Hydrogen Z=1, lowest priority)")
    print("\nStep 2: Determine configuration from 3D structure.")
    print("In the drawing, the lowest priority group (-H) is in the plane.")
    print("To solve this, we mentally swap the group in the back (-CH3, priority 3) with the group in the plane (-H, priority 4).")
    print("After the swap, -H is in the back. We trace the path from priority 1 -> 2 -> 3.")
    print("The path from -OH (1) -> -C(Et)... (2) -> -CH3 (3) is now Clockwise, which indicates (R).")
    print("However, because we performed one swap, the true configuration is the opposite.")
    c2_config = "S"
    print(f"Therefore, the configuration at C2 is (R) -> reversed -> {c2_config}\n")

    print("Analyzing the stereocenter at C3 (center, with ethyl group):")
    print("Step 1: Assign CIP priorities.")
    print("1: -C2(OH)... (the carbon on the left is bonded to an Oxygen)")
    print("2: -C4(Me)... (the carbon on the right is bonded to C, C, H)")
    print("3: -CH2CH3 (the ethyl group's carbon is bonded to C, H, H)")
    print("4: -H (lowest priority)")
    print("\nStep 2: Determine configuration from 3D structure.")
    print("The lowest priority group (-H) is pointing to the back (dashed wedge).")
    print("We can directly trace the path from priority 1 -> 2 -> 3.")
    print("The path from -C2(OH)... (1) -> -C4(Me)... (2) -> -CH2CH3 (3) is Clockwise.")
    c3_config = "R"
    print(f"Therefore, the configuration at C3 is {c3_config}\n")

    print("Analyzing the stereocenter at C4 (rightmost):")
    print("Step 1: Assign CIP priorities.")
    print("1: -CH2NH2 (the carbon is bonded to a Nitrogen, Z=7)")
    print("2: -C3(Et)... (the carbon is bonded to two other carbons, Z=6)")
    print("3: -CH3 (the carbon is bonded to three hydrogens, Z=1)")
    print("4: -H (lowest priority)")
    print("\nStep 2: Determine configuration from 3D structure.")
    print("The drawing places the lowest priority group (-H) in the front (solid wedge).")
    print("We trace the path from priority 1 -> 2 -> 3.")
    print("The path from -CH2NH2 (1) -> -C3(Et)... (2) -> -CH3 (3) is Counter-clockwise, which indicates (S).")
    print("Because the lowest priority group is in the front, we must reverse the result.")
    c4_config = "R"
    print(f"Therefore, the configuration at C4 is (S) -> reversed -> {c4_config}\n")
    
    print("--- Final Result ---")
    print(f"The absolute configuration of the molecule is ({2}{c2_config}, {3}{c3_config}, {4}{c4_config}).")

solve_stereochemistry()