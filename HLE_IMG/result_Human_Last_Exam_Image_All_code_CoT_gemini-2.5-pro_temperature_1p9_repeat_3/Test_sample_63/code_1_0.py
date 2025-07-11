import collections

def analyze_nmr_spectrum():
    """
    Analyzes an NMR spectrum and identifies the corresponding chemical structure.
    """
    # Step 1: Analyze the observed 1H NMR spectrum from the image.
    # The dictionary keys are the approximate chemical shifts (ppm).
    # The values are tuples of (Integration, Multiplicity).
    observed_spectrum = {
        1.1: ('6H', 'triplet'),
        2.3: ('3H', 'singlet'),
        2.7: ('4H', 'quartet'),
        3.4: ('2H', 'singlet'),
        7.1: ('3H', 'multiplet'), # Approximating integral for 3 aromatic protons
        9.0: ('1H', 'broad singlet')
    }

    print("--- Analysis of the Observed 1H NMR Spectrum ---")
    print("A detailed analysis of the signals reveals the following features:")
    print("- A quartet at ~2.7 ppm (4H) and a triplet at ~1.1 ppm (6H) strongly indicate the presence of two equivalent ethyl groups attached to a nitrogen: -N(CH2CH3)2.")
    print("- A broad singlet far downfield at ~9.0 ppm (1H) is characteristic of an amide proton (-NH-C=O).")
    print("- A singlet at ~3.4 ppm (2H) corresponds to an isolated methylene group (-CH2-) with no adjacent protons.")
    print("- A singlet at ~2.3 ppm (3H) corresponds to an isolated methyl group (-CH3).")
    print("- A multiplet around ~7.1 ppm corresponds to protons on an aromatic ring.\n")

    # Step 2 & 3: Evaluate each candidate structure.
    print("--- Evaluating Candidate Structures ---")

    # Candidate A-G and B-G
    print("Analysis of A-G and B-G:")
    print("Both structures contain a dimethylamino group, -N(CH3)2. This group would produce a single 6H singlet.")
    print("The observed spectrum does not have a 6H singlet; instead, it shows signals for two ethyl groups.")
    print("Result: A-G and B-G are incorrect.\n")

    # Candidate D-L
    print("Analysis of D-L:")
    print("This structure is 2-(diethylamino)-N-(2,4-dimethylphenyl)acetamide or similar (has two methyls on the ring).")
    print("It has two methyl groups on the aromatic ring, which would integrate to 6H in total (either as one 6H singlet if equivalent, or two 3H singlets if not).")
    print("The observed spectrum shows only one methyl singlet at ~2.3 ppm with an integration corresponding to 3H.")
    print("Result: D-L is incorrect.\n")
    
    # Candidate C-L
    print("Analysis of C-L:")
    print("This structure is 2-(diethylamino)-N-(2-methylphenyl)acetamide.")
    print("Let's compare its predicted signals to the observed spectrum:")
    
    # The 'equation' format: assigning each observed peak to a part of molecule C-L.
    print("\n--- Final Assignment for Structure C-L ---")
    print("Structure C-L matches the observed spectrum perfectly. Here is the peak-by-peak assignment:")
    print(f"Observed peak at 9.0 ppm = 1H broad singlet = Amide proton (-NH-C=O)")
    print(f"Observed peak at 7.1 ppm = 3H multiplet   = Protons on the substituted benzene ring")
    print(f"Observed peak at 3.4 ppm = 2H singlet     = Methylene protons in -CO-CH2-N=")
    print(f"Observed peak at 2.7 ppm = 4H quartet     = Methylene protons in -N-(CH2-CH3)2")
    print(f"Observed peak at 2.3 ppm = 3H singlet     = Methyl protons of Ar-CH3 on the ring")
    print(f"Observed peak at 1.1 ppm = 6H triplet     = Methyl protons in -N-(CH2-CH3)2")
    
    print("\nResult: C-L is the correct structure.\n")


if __name__ == '__main__':
    analyze_nmr_spectrum()
    print("<<<C>>>")
