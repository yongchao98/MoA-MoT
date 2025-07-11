def analyze_nmr_spectrum():
    """
    Analyzes an NMR spectrum to identify the most likely molecular structure
    from a list of candidates by following a logical deduction process.
    """
    print("### NMR Spectrum Analysis ###\n")

    # Define key features from the observed spectrum
    observed_triplet_ppm = 1.1
    observed_quartet_ppm = 2.5
    observed_aromatic_methyl_singlet_ppm = 2.3
    
    # Define features of candidate molecules
    candidates = {
        'A-G': {'has_ethyl_group': False},
        'B-G': {'has_ethyl_group': False},
        'C-L': {'has_ethyl_group': True, 'aromatic_methyl_protons': 3},
        'D-L': {'has_ethyl_group': True, 'aromatic_methyl_protons': 6}
    }
    
    # --- Step 1 & 2: Check for Ethyl Group and First Elimination ---
    print("Step 1: Identify key patterns in the spectrum.")
    print(f"The spectrum displays a triplet signal at ~{observed_triplet_ppm} ppm and a quartet signal at ~{observed_quartet_ppm} ppm.")
    print("This pattern is the classic signature of an ethyl group (-CH2-CH3).\n")
    
    possible_candidates = []
    print("Step 2: Eliminate candidates lacking an ethyl group.")
    for name, features in candidates.items():
        if features['has_ethyl_group']:
            print(f"- Candidate {name} contains diethylamino groups and is a possible match.")
            possible_candidates.append(name)
        else:
            print(f"- Candidate {name} does not contain ethyl groups and is eliminated.")
            
    print(f"\nThis narrows down the possibilities to: {', '.join(possible_candidates)}.\n")

    # --- Step 3 & 4: Differentiate remaining candidates using integration ---
    print("Step 3: Compare the remaining structures, C-L and D-L.")
    protons_c = candidates['C-L']['aromatic_methyl_protons']
    protons_d = candidates['D-L']['aromatic_methyl_protons']
    print(f"The primary difference is the number of protons from methyl groups on the aromatic ring.")
    print(f"- Structure C-L has one aromatic methyl group ({protons_c}H).")
    print(f"- Structure D-L has two aromatic methyl groups ({protons_d}H).\n")
    
    print("Step 4: Use integration ratios to identify the correct structure.")
    # In both C-L and D-L, the N-diethyl group gives a triplet for the -CH3 part.
    ethyl_methyl_protons = 6  # from -N(CH2-CH3)2
    print(f"We can use the triplet at ~{observed_triplet_ppm} ppm as a reference. It represents the two methyls of the ethyl groups, integrating to {ethyl_methyl_protons}H.")
    print(f"The singlet for the aromatic methyl group(s) is observed at ~{observed_aromatic_methyl_singlet_ppm} ppm.")

    print("\nLet's test the two hypotheses:")
    
    # Hypothesis for C-L
    print("  - If the structure is C-L:")
    print(f"    The ratio of protons for (aromatic CH3) to (ethyl CH3) should be {protons_c}H : {ethyl_methyl_protons}H.")
    print(f"    The final equation for the integration ratio is: {protons_c} : {ethyl_methyl_protons} = 1 : 2")

    # Hypothesis for D-L
    print("  - If the structure is D-L:")
    print(f"    The ratio of protons for (aromatic CH3) to (ethyl CH3) should be {protons_d}H : {ethyl_methyl_protons}H.")
    print(f"    The final equation for the integration ratio is: {protons_d} : {ethyl_methyl_protons} = 1 : 1")
    
    # --- Step 5: Final Conclusion ---
    print("\nStep 5: Conclusion from visual inspection.")
    print(f"Observing the spectrum, the area of the singlet at {observed_aromatic_methyl_singlet_ppm} ppm is clearly smaller than the area of the triplet at {observed_triplet_ppm} ppm.")
    print("The 1:1 integration ratio for D-L is incorrect. The 1:2 ratio fits the visual data well.")
    print("\nTherefore, the evidence strongly indicates that the molecule is C-L.")


# Execute the analysis function
analyze_nmr_spectrum()