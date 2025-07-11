def determine_configuration():
    """
    This function outlines the step-by-step determination of the absolute configuration
    for the three chiral centers in the provided molecule.
    """

    print("Step 1: Identify the chiral centers and their substituents.")
    print("The molecule has three chiral centers, labeled C2, C3, and C4 from left to right.")
    print(" - C2 is bonded to: -OH, -CH3, -H, and the rest of the molecule (C3).")
    print(" - C3 is bonded to: -CH2CH3 (ethyl), -H, and the carbons C2 and C4.")
    print(" - C4 is bonded to: -CH2NH2, -CH3, -H, and C3.")
    print("-" * 30)

    print("Step 2: Determine the configuration of C4 (the rightmost chiral center).")
    print("Substituents on C4: -CH2NH2, -C3_group, -CH3, -H.")
    print("Priorities (CIP Rules):")
    print("1: -CH2NH2 (due to N atom)")
    print("2: -C3_group (complex carbon group)")
    print("3: -CH3")
    print("4: -H")
    print("Spatial arrangement: -CH2NH2 is wedge (front), -CH3 is dash (back). The -H is in the plane.")
    print("To determine configuration, we can rotate the molecule to place H (priority 4) in the back. A 90-degree rotation around the C3-C4 bond axis places H in the back, -CH2NH2 pointing up, and -CH3 pointing down.")
    print("The sequence 1 -> 2 -> 3 (-CH2NH2 -> -C3_group -> -CH3) is counter-clockwise.")
    c4_config = "S"
    print(f"Result for C4: {c4_config}")
    print("-" * 30)

    print("Step 3: Determine the configuration of C3 (the middle chiral center).")
    print("Substituents on C3: -C2_group, -C4_group, -CH2CH3, -H.")
    print("Priorities (CIP Rules):")
    print("1: -C2_group (bonded to O atom)")
    print("2: -C4_group (bonded to a C which is bonded to N)")
    print("3: -CH2CH3 (ethyl)")
    print("4: -H")
    print("Spatial arrangement: -H is implied to be on a dash (back).")
    print("With the lowest priority group in the back, we trace the path 1 -> 2 -> 3.")
    print("The sequence (-C2_group -> -C4_group -> -CH2CH3) is clockwise.")
    c3_config = "R"
    print(f"Result for C3: {c3_config}")
    print("-" * 30)
    
    print("Step 4: Determine the configuration of C2 (the leftmost chiral center).")
    print("The drawing for C2 is ambiguous, showing two dashed bonds, which is impossible.")
    print("We assume a typo: the -OH group should be a wedge and the -CH3 group a dash.")
    print("Substituents on C2: -OH, -C3_group, -CH3, -H.")
    print("Priorities (CIP Rules):")
    print("1: -OH")
    print("2: -C3_group")
    print("3: -CH3")
    print("4: -H")
    print("Assumed spatial arrangement: -OH is wedge (front), -CH3 is dash (back). -H is in the plane.")
    print("With P1 (wedge) and P3 (dash), we can see that the sequence 1 -> 2 -> 3 is clockwise.")
    c2_config = "R"
    print(f"Result for C2: {c2_config}")
    print("-" * 30)

    print("Final Conclusion:")
    print(f"The absolute configuration of the molecule is {c2_config} at C2, {c3_config} at C3, and {c4_config} at C4.")
    print(f"Final Answer: 2{c2_config}, 3{c3_config}, 4{c4_config}")


if __name__ == '__main__':
    determine_configuration()