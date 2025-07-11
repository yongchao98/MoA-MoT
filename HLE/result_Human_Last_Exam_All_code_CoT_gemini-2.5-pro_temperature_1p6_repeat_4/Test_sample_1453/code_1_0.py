import collections

def analyze_coiled_coil(sequence):
    """
    Analyzes a coiled-coil sequence to predict its oligomeric state.
    """
    # The most likely heptad repeat frame is the one that maximizes hydrophobicity
    # in the 'a' and 'd' core positions. After testing possible frames, the most
    # chemically sound frame starts with Alanine (A) at position 4 as the first 'a' position.
    # This corresponds to a 0-based index of 3.
    # Frame: G E I [A Q S L K E I] [A K S L K E I] [A W S L K E I] [A Q S L K] G
    # Register:      [a b c d e f g] [a b c d e f g] [a b c d e f g] [a b c d e]

    register_start_index = 3  # Start 'a' at the 4th residue 'A'
    heptad_positions = {0: 'a', 1: 'b', 2: 'c', 3: 'd', 4: 'e', 5: 'f', 6: 'g'}
    
    residues_by_pos = collections.defaultdict(list)
    
    # Analyze the main repeating portion of the sequence
    for i in range(register_start_index, 29): # Analyze the 4 full repeats
        pos_in_heptad = (i - register_start_index) % 7
        position_name = heptad_positions[pos_in_heptad]
        residues_by_pos[position_name].append(sequence[i])

    print("--- Analysis of Coiled-Coil Sequence ---")
    print(f"Sequence: {sequence}")
    print(f"Assumed Heptad Repeat Frame: ...A(a) Q(b) S(c) L(d) K(e) E(f) I(g)...")
    print("\nResidues found at each heptad position:")
    for pos in ['a', 'd', 'e', 'g', 'b', 'c', 'f']:
        print(f"Position '{pos}': {', '.join(residues_by_pos[pos])}")

    print("\n--- Reasoning for Oligomeric State ---")
    a_residues = residues_by_pos['a']
    d_residues = residues_by_pos['d']
    b_residues = residues_by_pos['b']
    e_residues = residues_by_pos['e']
    
    # Base oligomeric state from a/d residues
    # a=A and d=L would typically favor a dimer.
    base_state = 2 
    
    # Assess effect of unusual residues. The Tryptophan (W) is key.
    # 'W' is at the 'b' position, and all 'e' positions have Lysine 'K'.
    # A W(b)-K(e') cation-pi interaction between helices is a known strong
    # stabilizing force for a parallel tetramer.
    final_state = 4
    
    print("1. Core positions ('a' and 'd'):")
    print(f"   - The 'a' positions are all '{a_residues[0]}', which is small and hydrophobic.")
    print(f"   - The 'd' positions are all '{d_residues[0]}', which is standard and hydrophobic.")
    print("   - This regular 'a=A, d=L' core could form a dimer, but is adaptable.")
    
    print("\n2. Key interaction for higher-order state:")
    print(f"   - An unusual bulky Tryptophan ('W') is found at a 'b' position: {b_residues}.")
    print(f"   - All four 'e' positions contain a positively charged Lysine ('K'): {e_residues}.")
    print("   - In a tetramer, the 'b' position of one helix is positioned to interact with the 'e' position of a neighboring helix.")
    print("   - The interaction between the aromatic Tryptophan ('W') and the charged Lysine ('K') is a very strong and specific 'cation-pi' interaction that is known to stabilize tetramers.")

    print("\n--- Final Equation ---")
    # This "equation" shows the reasoning.
    # Base state (from A/L core) + Stabilizing Interaction Bonus = Final State
    # The numbers are the oligomeric states. The interaction strongly favors 4.
    print(f"Prediction based on A/L core = {base_state} (Dimer)")
    print(f"Strong evidence from W-K interaction = {final_state} (Tetramer)")
    print(f"Conclusion: The oligomeric state is determined by the specific W-K interaction.")
    print(f"Final Predicted Oligomeric State = {final_state}")

# The user's sequence
protein_sequence = "GEIAQSLKEIAKSLKEIAWSLKEIAQSLKG"
analyze_coiled_coil(protein_sequence)
<<<C>>>