def determine_configuration():
    """
    Determines and explains the absolute configuration of the molecule shown in the image.
    The molecule is identified as 5-amino-3-ethyl-4-methylpentan-2-ol.
    """
    print("Determining the absolute configuration of 5-amino-3-ethyl-4-methylpentan-2-ol.")
    print("The carbon chain is numbered to give the -OH group the lowest number (C2).\n")

    # --- Analysis of Chiral Center C2 ---
    print("--- Analysis of C2 ---")
    print("The drawing at C2 is ambiguous. It's assumed to represent an anti-conformation relative to C3.")
    print("This means: -OH is back (dashed) and -CH3 is forward (wedged). The H atom is in the plane.")
    print("Priority of groups attached to C2 (Cahn-Ingold-Prelog):")
    c2_priorities = {"-OH": 1, "-C3(group)": 2, "-CH3": 3, "-H": 4}
    print(f"  1: {list(c2_priorities.keys())[0]}\n  2: {list(c2_priorities.keys())[1]}\n  3: {list(c2_priorities.keys())[2]}\n  4: {list(c2_priorities.keys())[3]}")
    # Using the swap method: H(4) is in the plane. Swap with -OH(1) which is back.
    # H is now back. The new arrangement of 1,2,3 is CCW -> S.
    # Original configuration is the opposite -> R.
    c2_config = "R"
    print(f"\nAfter analysis (e.g., using the swap method), the configuration is determined to be clockwise.")
    print(f"Configuration at C2 is {c2_config}.\n")

    # --- Analysis of Chiral Center C3 ---
    print("--- Analysis of C3 ---")
    print("Priority of groups attached to C3 (Cahn-Ingold-Prelog):")
    c3_priorities = {"-C2(group)": 1, "-C4(group)": 2, "-Ethyl": 3, "-H": 4}
    print(f"  1: {list(c3_priorities.keys())[0]}\n  2: {list(c3_priorities.keys())[1]}\n  3: {list(c3_priorities.keys())[2]}\n  4: {list(c3_priorities.keys())[3]}")
    print("\nThe lowest priority group (-H) is in the back (dashed).")
    print("The path from priority 1 -> 2 -> 3 is clockwise.")
    c3_config = "R"
    print(f"Configuration at C3 is {c3_config}.\n")

    # --- Analysis of Chiral Center C4 ---
    print("--- Analysis of C4 ---")
    print("Priority of groups attached to C4 (Cahn-Ingold-Prelog):")
    # C5 is the -CH2NH2 group
    c4_priorities = {"-C5(CH2NH2)": 1, "-C3(group)": 2, "-CH3": 3, "-H": 4}
    print(f"  1: {list(c4_priorities.keys())[0]}\n  2: {list(c4_priorities.keys())[1]}\n  3: {list(c4_priorities.keys())[2]}\n  4: {list(c4_priorities.keys())[3]}")
    print("\nThe lowest priority group (-H) is in the front (wedged).")
    print("The path from priority 1 -> 2 -> 3 is counter-clockwise, which suggests S.")
    print("However, since the lowest priority group is in the front, the assignment is reversed.")
    c4_config = "R"
    print(f"Configuration at C4 is {c4_config}.\n")

    # --- Final Result ---
    print("--- Final Result ---")
    print(f"The complete absolute configuration of the molecule is:")
    # The final equation is built here, showing each number and letter
    final_config_str = f"({2}{c2_config}, {3}{c3_config}, {4}{c4_config})"
    print(final_config_str)


if __name__ == '__main__':
    determine_configuration()