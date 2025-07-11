def analyze_nmr_spectrum():
    """
    Analyzes an H-NMR spectrum to determine the most likely molecular structure
    from a list of candidates.
    """
    print("Step 1: Analyzing the provided 1H NMR Spectrum")
    print("-------------------------------------------------")
    print("We identify the following signals:")
    print(" - Signal at ~9.0 ppm: Broad singlet, relative integral = 1H. Characteristic of an amide N-H proton.")
    print(" - Signal at ~7.1 ppm: Multiplet, relative integral = 3H. Characteristic of a tri-substituted benzene ring.")
    print(" - Signal at ~3.4 ppm: Singlet, relative integral = 2H. Characteristic of a CH2 group with no adjacent protons, likely between a C=O and N atom.")
    print(" - Signal at ~2.8 ppm: Quartet, relative integral = 4H. Characteristic of two equivalent CH2 groups adjacent to CH3 groups.")
    print(" - Signal at ~2.3 ppm: Singlet, relative integral = 6H. Characteristic of two equivalent CH3 groups with no adjacent protons.")
    print(" - Signal at ~1.1 ppm: Triplet, relative integral = 6H. Characteristic of two equivalent CH3 groups adjacent to CH2 groups.")

    print("\nStep 2: Evaluating the Candidate Structures")
    print("-------------------------------------------")
    print(" - Candidates A-G and B-G (Gramine/DMT): These are indole derivatives without an amide bond or ethyl groups. They are inconsistent with the observed quartet-triplet pattern and amide N-H signal. They are eliminated.")
    print(" - Candidates C-L and D-L: Both are N-(dimethylphenyl)-2-(diethylamino)acetamide derivatives. They both have the necessary fragments (amide, diethylamino group, phenyl ring) to match most signals.")

    print("\nStep 3: Differentiating between C-L and D-L")
    print("----------------------------------------------")
    print("The key difference is the aromatic substitution pattern.")
    print(" - Structure D-L is N-(2,5-dimethylphenyl)... The methyl groups at positions 2 and 5 are not chemically equivalent and should produce TWO separate 3H singlets.")
    print(" - Structure C-L is N-(2,6-dimethylphenyl)... (Lidocaine). Due to molecular symmetry, the methyl groups at positions 2 and 6 are chemically EQUIVALENT. They should produce ONE 6H singlet.")

    print("\nStep 4: Final Conclusion")
    print("-------------------------")
    print("The spectrum shows ONE singlet at ~2.3 ppm with an integration corresponding to 6H. This is only consistent with the symmetrical structure of C-L.")
    print("\nFinal match of spectrum to Structure C-L (Lidocaine):")
    print(f" - ~9.0 ppm (s, 1H)  ->  Amide N-H")
    print(f" - ~7.1 ppm (m, 3H)  ->  Aromatic H's")
    print(f" - ~3.4 ppm (s, 2H)  ->  -CO-CH2-N-")
    print(f" - ~2.3 ppm (s, 6H)  ->  Two equivalent Ar-CH3 groups")
    print(f" - ~2.8 ppm (q, 4H) and ~1.1 ppm (t, 6H) -> Diethylamino group -N(CH2CH3)2")
    print("\nThe spectrum perfectly matches structure C-L.")

# Run the analysis
analyze_nmr_spectrum()