def analyze_nmr_spectrum():
    """
    Analyzes an H-NMR spectrum to identify the correct molecular structure from a list of candidates.
    The function will print a step-by-step analysis and the final conclusion.
    """

    # Step 1: Analyze the observed spectrum
    print("Step 1: Analysis of the Provided 1H NMR Spectrum")
    print("-------------------------------------------------")
    print("The spectrum shows the following signals:")
    print("- A triplet at ~1.2 ppm (let's denote its integration as 6H)")
    print("- A singlet at ~2.3 ppm (integration ~3H)")
    print("- A quartet at ~2.8 ppm (integration ~4H)")
    print("- A singlet at ~3.3 ppm (integration ~2H)")
    print("- A multiplet in the aromatic region at ~7.1-7.2 ppm (integration ~4H)")
    print("- A broad singlet far downfield at ~9.0 ppm (integration ~1H)")
    print("\nThe quartet (~2.8 ppm) and triplet (~1.2 ppm) are characteristic of two ethyl groups (-CH2CH3).\n")

    # Step 2: Evaluate each candidate structure
    print("Step 2: Evaluating Candidate Structures")
    print("---------------------------------------")

    # Candidate A-G and B-G
    print("Candidates A-G and B-G:")
    print("  - These structures contain a dimethylamino group [-N(CH3)2], not diethylamino groups.")
    print("  - They lack the ethyl groups needed to produce the observed quartet and triplet.")
    print("  - Conclusion: A-G and B-G are incorrect.\n")

    # Candidate D-L (N-(2,6-dimethylphenyl)...)
    print("Candidate D-L:")
    print("  - This structure has two methyl groups on the aromatic ring (2,6-dimethyl).")
    print("  - This would produce a singlet with an integration of 6H for these two equivalent methyls.")
    print("  - The spectrum shows a singlet at ~2.3 ppm with an integration of only 3H.")
    print("  - Also, it only has 3 aromatic protons, while the spectrum suggests 4.")
    print("  - Conclusion: D-L is incorrect.\n")
    
    # Candidate C-L (N-(2-methylphenyl)...)
    print("Candidate C-L:")
    print("  - Let's predict the signals for this structure:")
    print("    - Amide N-H: 1 proton, gives a broad singlet at ~9.0 ppm. -> MATCHES (1H)")
    print("    - Aromatic Ring: 4 non-equivalent protons, giving a multiplet at ~7.1-7.2 ppm. -> MATCHES (4H)")
    print("    - Methylene (-CO-CH2-N): 2 protons, no adjacent protons, gives a singlet at ~3.3 ppm. -> MATCHES (2H)")
    print("    - Diethylamino Group [-N(CH2CH3)2]:")
    print("      - Methylene part (-CH2-): 4 protons total, adjacent to 3 protons, gives a quartet at ~2.8 ppm. -> MATCHES (4H)")
    print("      - Methyl part (-CH3): 6 protons total, adjacent to 2 protons, gives a triplet at ~1.2 ppm. -> MATCHES (6H)")
    print("    - Aromatic Methyl (Ar-CH3): 3 protons, no adjacent protons, gives a singlet at ~2.3 ppm. -> MATCHES (3H)")
    print("  - Conclusion: All signals in the spectrum are perfectly explained by structure C-L.\n")

    # Final Conclusion
    print("Step 3: Final Conclusion")
    print("------------------------")
    print("Based on the analysis, the spectrum corresponds to structure C-L, which is 2-(diethylamino)-N-(2-methylphenyl)acetamide.")
    print("Final Equation (Proton Count): 1 (NH) + 4 (Ar-H) + 2 (-COCH2-) + 4 (-NCH2-) + 3 (Ar-CH3) + 6 (-CH2CH3) = 20 total protons.")


# Run the analysis
analyze_nmr_spectrum()
