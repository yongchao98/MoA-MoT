def identify_hydrocarbon():
    """
    Identifies a hydrocarbon from its molecular formula and 13C NMR spectrum.
    """
    # Provided Data
    molecular_formula = "C7H14"
    nmr_data = {
        145: 's (singlet, quaternary C)',
        112: 't (triplet, CH2)',
        48: 't (triplet, CH2)',
        27: 'd (doublet, CH)',
        22: 'q (quartet, CH3)',
        21: 'q (quartet, CH3)'
    }

    print("### Step-by-Step Analysis ###\n")

    # Step 1: Analyze the Molecular Formula
    print("1. Degree of Unsaturation (DBE):")
    carbons = 7
    hydrogens = 14
    dbe = int(carbons - (hydrogens / 2) + 1)
    print(f"   - The formula {molecular_formula} has a DBE of {dbe}.")
    print("   - This indicates the presence of one double bond or one ring.\n")

    # Step 2: Analyze the 13C NMR Spectrum
    print("2. Analysis of ¹³C NMR Signals:")
    print("   - Signals at 145(s) and 112(t) are in the alkene region (100-160 ppm), which confirms the double bond.")
    print("   - The signal at 145(s) corresponds to a quaternary carbon (>C=) in the double bond.")
    print("   - The signal at 112(t) corresponds to a terminal methylene carbon (=CH₂) in the double bond.")
    print("   - These two signals strongly indicate a >C=CH₂ fragment.\n")
    print("   - The spectrum shows 6 signals for a 7-carbon molecule. This means two carbons must be chemically equivalent.")
    print("   - A count of hydrogens from the multiplicities (assuming one 'q' signal is for two carbons) is: 0 (from 145s) + 2 (from 112t) + 2 (from 48t) + 1 (from 27d) + 3 (from one q) + 2*3 (from the other q) = 14 H. This matches the formula, confirming two equivalent CH₃ groups.\n")

    # Step 3: Identify and Assemble Fragments
    print("3. Assembling the Structure:")
    print("   Based on the signals, we can identify three structural units:")
    print("   - Fragment A (from 145s, 112t): >C=CH₂")
    print("   - Fragment B (from 48t, 22q): An ethyl group (-CH₂CH₃)")
    print("   - Fragment C (from 27d, 21q x2): An isopropyl group (-CH(CH₃)₂)\n")
    print("   The quaternary carbon from Fragment A must connect to the ethyl and isopropyl groups.")
    print("   Structure: (CH₃)₂CH - C(CH₂CH₃) = CH₂\n")

    # Step 4: Determine the IUPAC Name
    print("4. IUPAC Nomenclature:")
    print("   - To name the structure, we must find the longest carbon chain containing the double bond.")
    print("   - Choosing the chain that goes through the isopropyl group gives a 4-carbon parent chain (methylbut-1-ene) with two substituents (ethyl at C2, methyl at C3).")
    print("   - Choosing the chain that goes through the ethyl group gives a 4-carbon parent chain (but-1-ene) with one substituent (isopropyl at C2).")
    print("   - IUPAC rules state that for chains of equal length, the one with the most substituents is chosen.")
    print("   - Therefore, the correct parent chain is 3-methylbut-1-ene, with an ethyl group at position 2.")
    final_name = "2-ethyl-3-methylbut-1-ene"
    print(f"\nFinal IUPAC Name Equation: {final_name}")

identify_hydrocarbon()
<<<2-ethyl-3-methylbut-1-ene>>>