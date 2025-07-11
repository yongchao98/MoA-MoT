import sys

def predict_oligomeric_state(sequence):
    """
    Predicts the oligomeric state of a coiled-coil sequence.
    """
    # Kyte-Doolittle hydrophobicity scale for scoring
    hydrophobicity = { 'A': 1.8, 'C': 2.5, 'D': -3.5, 'E': -3.5, 'F': 2.8,
                       'G': -0.4, 'H': -3.2, 'I': 4.5, 'K': -3.9, 'L': 3.8,
                       'M': 1.9, 'N': -3.5, 'P': -1.6, 'Q': -3.5, 'R': -4.5,
                       'S': -0.8, 'T': -0.7, 'V': 4.2, 'W': -0.9, 'Y': -1.3 }

    print(f"Analyzing sequence: {sequence}\n")
    print("Step 1: Identifying the most likely heptad repeat ('abcdefg') alignment.")
    print("This is done by finding the frame that maximizes the hydrophobicity of the core 'a' and 'd' positions.\n")

    best_offset = -1
    max_score = -sys.maxsize
    best_a = []
    best_d = []
    
    # Test all 7 possible offsets for the heptad repeat
    for offset in range(7):
        a_residues = []
        d_residues = []
        current_score = 0.0
        
        for i, aa in enumerate(sequence):
            # register 'a' is index 0, 'b' is 1, ..., 'g' is 6
            register_index = (i + offset) % 7
            
            if register_index == 0:  # 'a' position
                a_residues.append(aa)
                current_score += hydrophobicity.get(aa, 0)
            elif register_index == 3:  # 'd' position
                d_residues.append(aa)
                current_score += hydrophobicity.get(aa, 0)
        
        if current_score > max_score:
            max_score = current_score
            best_offset = offset
            best_a = a_residues
            best_d = d_residues
            
    heptad_map = ['a', 'b', 'c', 'd', 'e', 'f', 'g']
    start_register = heptad_map[(0 + best_offset) % 7]
    
    print("Step 2: Best alignment found.")
    print(f"The best alignment starts with the first residue ('G') in the '{start_register}' position.")
    print(f"This alignment yields a maximum core hydrophobicity score of {max_score:.1f}.")
    print("The identified core residues are:")
    print(f"  'a' positions: {best_a}")
    print(f"  'd' positions: {best_d}\n")
    
    # Step 3: Analyze the pattern and predict
    print("Step 3: Analyzing core composition to predict oligomeric state.")
    oligomeric_state = "Unknown"
    reasoning = "The pattern does not match common motifs."

    # Convert to sets for easier comparison, ignoring potential outliers at ends
    a_set = set(best_a[:4])
    d_set = set(best_d[:4])
    
    if a_set == {'L'} and d_set == {'I'}:
        oligomeric_state = 2
        reasoning = ("This pattern, with Leucine (L) at the 'a' positions and Isoleucine (I) at the 'd' positions, "
                     "promotes a dimeric coiled-coil. The bulky, beta-branched side chain of Isoleucine at the 'd' "
                     "position packs efficiently against the Leucine at the 'a' position of the partner helix in a dimer.")
    elif a_set == {'I'} and d_set == {'L'}:
        oligomeric_state = 3
        reasoning = ("This pattern, with Isoleucine (I) at the 'a' positions and Leucine (L) at the 'd' positions, "
                     "is a classic hallmark of a trimeric coiled-coil. The bulky side chain of Isoleucine fits best in the "
                     "more spacious core created by three helices at the 'a' positions.")
    elif a_set == {'A'} and d_set == {'L'}:
         oligomeric_state = 2
         reasoning = ("This is a canonical 'leucine zipper' motif. The small Alanine (A) at 'a' and Leucine (L) at 'd' "
                      "pack together very efficiently to form a stable dimer.")

    print(f"Reasoning: {reasoning}")
    print("\n--- Conclusion ---")
    print(f"The predicted oligomeric state is: {oligomeric_state}")

# The sequence provided by the user
protein_sequence = "GEIAQSLKEIAKSLKEIAWSLKEIAQSLKG"
predict_oligomeric_state(protein_sequence)
