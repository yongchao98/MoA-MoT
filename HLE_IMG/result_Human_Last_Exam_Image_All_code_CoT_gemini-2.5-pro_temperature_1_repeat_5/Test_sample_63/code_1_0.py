def analyze_nmr_spectrum():
    """
    Analyzes the provided 1H NMR spectrum and matches it to the correct chemical structure.
    """
    print("Step 1: Analyzing the signals from the 1H NMR spectrum.")
    print("---------------------------------------------------------")
    print("Signal at ~9.0 ppm: A broad singlet corresponding to 1 proton. This is characteristic of an amide N-H proton.")
    print("Signal at ~7.1 ppm: A complex multiplet corresponding to 4 protons. This is the aromatic region for a mono-substituted phenyl ring.")
    print("Signal at ~3.5 ppm: A sharp singlet corresponding to 2 protons. This indicates a CH2 group with no adjacent protons, likely -CO-CH2-N-.")
    print("Signal at ~2.7 ppm: A quartet corresponding to 4 protons. This indicates two equivalent CH2 groups next to CH3 groups.")
    print("Signal at ~2.3 ppm: A sharp singlet corresponding to 3 protons. This indicates a CH3 group with no adjacent protons, likely on the aromatic ring.")
    print("Signal at ~1.1 ppm: A triplet corresponding to 6 protons. This indicates two equivalent CH3 groups next to CH2 groups.")
    print("\nStep 2: Matching the spectral data to the candidate molecules.")
    print("------------------------------------------------------------")
    print("Candidates A-G and B-G are eliminated because they lack the ethyl groups (-CH2CH3) that produce the characteristic triplet and quartet signals.")
    print("Candidate D-L is eliminated because it has two methyl groups on the aromatic ring, which would produce a singlet for 6 protons at ~2.3 ppm, not 3 protons as seen in the spectrum.")
    print("\nStep 3: Confirming the match with Candidate C-L.")
    print("---------------------------------------------------")
    print("Structure C-L has:")
    print("- An amide N-H proton: Matches the 1H broad singlet at 9.0 ppm.")
    print("- Four aromatic protons: Matches the 4H multiplet at 7.1 ppm.")
    print("- A -CO-CH2-N- group: Matches the 2H singlet at 3.5 ppm.")
    print("- Two ethyl groups on a nitrogen: The -CH2- parts match the 4H quartet at 2.7 ppm and the -CH3 parts match the 6H triplet at 1.1 ppm.")
    print("- One methyl group on the aromatic ring: Matches the 3H singlet at 2.3 ppm.")
    print("\nConclusion: The spectrum perfectly matches the structure of C-L.")

analyze_nmr_spectrum()