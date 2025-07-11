def solve_nmr():
    """
    Analyzes an NMR spectrum to identify the corresponding molecular structure.
    This function will print the reasoning step-by-step.
    """

    # --- Spectrum Analysis ---
    print("Step 1: Analyzing the 1H NMR Spectrum")
    print("----------------------------------------")
    signals = {
        "8.8 ppm": {"protons": 1, "multiplicity": "broad singlet", "assignment": "Amide N-H"},
        "7.2 ppm": {"protons": "3 or 4", "multiplicity": "multiplet", "assignment": "Aromatic Ar-H"},
        "3.4 ppm": {"protons": 2, "multiplicity": "singlet", "assignment": "C(=O)-CH2-N"},
        "2.8 ppm": {"protons": 4, "multiplicity": "quartet", "assignment": "-N-(CH2CH3)2"},
        "2.3 ppm": {"protons": 3, "multiplicity": "singlet", "assignment": "Aromatic Ar-CH3"},
        "1.2 ppm": {"protons": 6, "multiplicity": "triplet", "assignment": "-N-(CH2CH3)2"},
    }

    print("Observed signals and their interpretation:")
    for shift, data in signals.items():
        print(f"- Peak at {shift}: A {data['multiplicity']} corresponding to {data['protons']}H. Likely from a {data['assignment']} group.")
    print("\n")

    # --- Candidate Analysis ---
    print("Step 2: Evaluating Candidate Structures")
    print("---------------------------------------")
    print("Candidates A-G and B-G: These structures lack ethyl groups (-CH2CH3). The spectrum clearly shows a 4H quartet and a 6H triplet characteristic of two ethyl groups. So, A-G and B-G are incorrect.")
    print("\nCandidates C-L and D-L: Both have the required N,N-diethyl acetamide side chain, matching the signals at 8.8, 3.4, 2.8, and 1.2 ppm.")
    print("The key difference is the aromatic ring substitution:")
    
    # Analysis for D-L
    aromatic_ch3_protons_D = 6
    aromatic_H_protons_D = 3
    print(f"- Structure D-L: Has two methyl groups on the ring. It should have a singlet for {aromatic_ch3_protons_D}H (Aromatic-CH3) and signals for {aromatic_H_protons_D}H in the aromatic region.")

    # Analysis for C-L
    aromatic_ch3_protons_C = 3
    aromatic_H_protons_C = 4
    print(f"- Structure C-L: Has one methyl group on the ring. It should have a singlet for {aromatic_ch3_protons_C}H (Aromatic-CH3) and signals for {aromatic_H_protons_C}H in the aromatic region.")
    print("\n")
    
    # --- Conclusion ---
    print("Step 3: Conclusion")
    print("------------------")
    observed_ar_ch3_protons = signals["2.3 ppm"]["protons"]
    print(f"The spectrum shows a singlet at 2.3 ppm integrating to {observed_ar_ch3_protons}H.")
    print(f"This observation ({observed_ar_ch3_protons}H) matches the predicted signal for structure C-L ({aromatic_ch3_protons_C}H) but not for structure D-L ({aromatic_ch3_protons_D}H).")
    print("Therefore, the most possible structure is C-L.")


solve_nmr()