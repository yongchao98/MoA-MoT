import collections

def predict_coiled_coil_state(sequences):
    """
    Predicts the oligomeric state of coiled-coil sequences based on
    residues at the 'a' and 'd' core positions of the heptad repeat.
    """
    
    oligomeric_states = []
    
    print("Analyzing coiled-coil sequences...")
    print("-" * 30)

    for i, seq in enumerate(sequences):
        # Assume a 'b' register, where the first residue is 'b' and the second is 'a'.
        # 'a' positions are at indices 1, 8, 15, 22, ...
        # 'd' positions are at indices 4, 11, 18, 25, ...
        a_positions = [seq[j] for j in range(1, len(seq), 7)]
        d_positions = [seq[j] for j in range(4, len(seq), 7)]
        
        a_counts = collections.Counter(a_positions)
        d_counts = collections.Counter(d_positions)

        state = "Unknown"
        reason = ""

        # Rule 1: A polar residue (N or Q) in the core is a dominant dimer-specifying feature.
        # This applies to Sequence 2 (N at 'a') and Sequence 4 (Q at 'd').
        if 'N' in a_positions or 'Q' in d_positions:
            state = 2
            reason = "A polar residue (N or Q) in a core position strongly specifies a dimer."
        
        # Rule 2: A core of Isoleucine (I) at both 'a' and 'd' positions specifies a tetramer.
        # This applies to Sequence 3.
        elif a_counts.get('I', 0) >= 3 and d_counts.get('I', 0) >= 3:
            state = 4
            reason = "A core of Isoleucine at both 'a' and 'd' positions favors a tetramer due to sterics."
            
        # Rule 3: A core of Isoleucine (I) at 'a' and Threonine (T) at 'd' is a known trimer motif.
        # This applies to Sequence 5.
        elif a_counts.get('I', 0) >= 3 and d_counts.get('T', 0) >= 3:
            state = 3
            reason = "A core of Isoleucine at 'a' and Threonine at 'd' is a known trimer-forming motif."

        # Rule 4: A core of Isoleucine (I) at 'a' and Alanine (A) at 'd' favors a dimer.
        # This applies to Sequence 1.
        elif a_counts.get('I', 0) >= 3 and d_counts.get('A', 0) >= 3:
            state = 2
            reason = "A core of Isoleucine at 'a' and the small Alanine at 'd' typically forms a dimer."
        
        oligomeric_states.append(state)
        
        print(f"Sequence {i+1}: {seq}")
        print(f"Core 'a' residues: {a_positions}")
        print(f"Core 'd' residues: {d_positions}")
        print(f"Prediction: {state} ({reason})")
        print("-" * 30)

    print("Final predicted oligomeric states:")
    # The problem asks to output the numbers in the final equation
    # We will format it as a simple list for clarity.
    print(oligomeric_states)

# The provided sequences
protein_sequences = [
    "EIAQALKEIAKALKEIAWALKEIAQALK", # Expected: 2
    "EIAALKQEIAALKKENAALKQEIAALKQ", # Expected: 2
    "EIAAIKQEIAAIKKEIAAIKWEIAAIKQ", # Expected: 4
    "EIQKQLKEIQKQLKEIQWQLKEIQKQLK", # Expected: 2
    "EIAQTLKEIAKTLKEIAWTLKEIAQTLK", # Expected: 3
]

predict_coiled_coil_state(protein_sequences)