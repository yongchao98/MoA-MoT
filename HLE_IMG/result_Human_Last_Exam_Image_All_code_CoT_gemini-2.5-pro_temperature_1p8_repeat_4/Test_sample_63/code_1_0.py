def analyze_nmr_spectrum():
    """
    Analyzes an NMR spectrum and compares it with candidate molecules
    to determine the most likely structure.
    """
    print("Step 1: Analyzing the experimental 1H NMR Spectrum")
    print("-----------------------------------------------------")
    spectrum_peaks = {
        "Peak 1": {"shift": "~8.9 ppm", "splitting": "Broad Singlet", "integration": "1H", "notes": "Likely an amide or amine N-H"},
        "Peak 2": {"shift": "~7.1 ppm", "splitting": "Multiplet", "integration": "~3H", "notes": "Aromatic protons (Ar-H)"},
        "Peak 3": {"shift": "~3.4 ppm", "splitting": "Singlet", "integration": "~2H", "notes": "Protons on carbon next to N or O, no adjacent H's"},
        "Peak 4": {"shift": "~2.8 ppm", "splitting": "Quartet", "integration": "~4H", "notes": "A -CH2- group next to a -CH3 group (part of an ethyl group)"},
        "Peak 5": {"shift": "~2.4 ppm", "splitting": "Singlet", "integration": "~3H", "notes": "A -CH3 group with no adjacent H's, likely on an aromatic ring or N"},
        "Peak 6": {"shift": "~1.2 ppm", "splitting": "Triplet", "integration": "~6H", "notes": "A -CH3 group next to a -CH2- group (part of an ethyl group)"},
    }
    for peak, data in spectrum_peaks.items():
        print(f"- {peak}: Shift={data['shift']}, Splitting={data['splitting']}. Notes: {data['notes']}")
    print("\nKey takeaway: The presence of a quartet and a triplet strongly suggests an ethyl group (-CH2CH3).\n")

    print("Step 2: Predicting spectra for candidate molecules")
    print("--------------------------------------------------")

    # Analysis of A-G
    print("Analyzing Candidate A-G:")
    print("  - Contains two -CH3 groups on a Nitrogen, expecting a 6H singlet.")
    print("  - Contains a -CH2- group, likely a singlet.")
    print("  - No ethyl group (-CH2CH3) present.")
    print("  - Verdict: Does NOT match the spectrum (missing quartet and triplet).\n")

    # Analysis of B-G
    print("Analyzing Candidate B-G:")
    print("  - Contains two -CH3 groups on a Nitrogen, expecting a 6H singlet.")
    print("  - No ethyl group (-CH2CH3) present.")
    print("  - Verdict: Does NOT match the spectrum (missing quartet and triplet).\n")

    # Analysis of D-L
    print("Analyzing Candidate D-L:")
    print("  - Aromatic protons (2H) are equivalent, so we expect a singlet in the aromatic region (7-8 ppm).")
    print("  - The experimental spectrum shows a multiplet, not a singlet.")
    print("  - Contains two equivalent aromatic -CH3 groups, expecting a 6H singlet.")
    print("  - The experimental spectrum has a smaller singlet at ~2.4 ppm.")
    print("  - Verdict: Does NOT match the spectrum (aromatic splitting and methyl integration are wrong).\n")
    
    # Analysis of C-L
    print("Analyzing Candidate C-L (Lidocaine analog):")
    print("Let's assign the protons of C-L to the observed peaks:")
    print(f"  1. Amide N-H (1H): Predicts a broad singlet around 8-9 ppm. -> MATCHES the peak at {spectrum_peaks['Peak 1']['shift']}.")
    print(f"  2. Aromatic Protons (3H): Non-equivalent, predicts a multiplet around 7-8 ppm. -> MATCHES the peak at {spectrum_peaks['Peak 2']['shift']}.")
    print(f"  3. Methylene C(O)-CH2-N (2H): No adjacent protons, predicts a singlet. -> MATCHES the peak at {spectrum_peaks['Peak 3']['shift']}.")
    print(f"  4. Aromatic Methyl Ar-CH3 (3H): No adjacent protons, predicts a singlet. -> MATCHES the peak at {spectrum_peaks['Peak 5']['shift']}.")
    print("  5. Diethylamino N(CH2CH3)2 group (4H + 6H):")
    print(f"     - Two N-CH2- groups (total 4H) are adjacent to CH3 groups. Predicts a quartet. -> MATCHES the peak at {spectrum_peaks['Peak 4']['shift']}.")
    print(f"     - Two -CH3 groups (total 6H) are adjacent to CH2 groups. Predicts a triplet. -> MATCHES the peak at {spectrum_peaks['Peak 6']['shift']}.")
    print("  - Verdict: The predicted spectrum for C-L is an excellent match for all observed signals.\n")
    
    print("Step 3: Conclusion")
    print("--------------------")
    print("Based on the detailed analysis, structure C-L is the only candidate whose predicted 1H NMR spectrum corresponds to the experimental data.")

if __name__ == "__main__":
    analyze_nmr_spectrum()