import collections

def predict_oligomeric_states(sequences):
    """
    Predicts the oligomeric state of coiled-coil sequences based on core residue analysis.

    This prediction uses a fixed heptad register g(abcdefg) and a set of rules
    based on the identity of residues at the 'a' and 'd' positions.
    """
    predictions = []
    
    print("Analysis:")
    for seq in sequences:
        # Step 1 & 2: Assume g(abcdefg) frame and extract a and d residues
        a_residues = [seq[i] for i in [1, 8, 15, 22]]
        d_residues = [seq[i] for i in [4, 11, 18, 25]]
        
        a_counts = collections.Counter(a_residues)
        d_counts = collections.Counter(d_residues)
        
        prediction = 0
        explanation = ""

        # Step 3: Apply prediction rules in a specific order
        # Rule 1: Polar 'K' in the 'd' position without a bulky 'W' suggests a tetramer.
        if 'K' in d_counts and 'W' not in d_counts:
            prediction = 4
            explanation = f"Detected 'K' at a 'd' position without 'W', predicting Tetramer ({prediction})."
        # Rule 2: A core composed entirely of 'T' at the 'd' positions is unusual and can support a trimer.
        elif d_counts.get('T') == 4:
            prediction = 3
            explanation = f"Detected all 'd' positions are 'T', which can stabilize a Trimer ({prediction})."
        # Rule 3: A bulky 'W' in the core can be accommodated, but often defaults to the simplest stable state, a dimer.
        elif 'W' in d_counts:
            prediction = 2
            explanation = f"Detected bulky 'W' at a 'd' position. This often results in a Dimer ({prediction})."
        # Rule 4: The canonical Isoleucine('a')/Leucine('d') pairing is a strong signal for a dimer.
        elif a_counts.most_common(1)[0][0] == 'I' and d_counts.most_common(1)[0][0] == 'L':
            prediction = 2
            explanation = f"Detected canonical 'I' at 'a' and 'L' at 'd' pairing, predicting Dimer ({prediction})."
        # Rule 5: A strong preference for 'I' at the 'a' position generally favors dimers.
        elif a_counts.most_common(1)[0][0] == 'I':
            prediction = 2
            explanation = f"Strong preference for 'I' at 'a' positions favors a Dimer ({prediction})."
        else:
            prediction = "Unknown"
            explanation = "Could not determine state with given rules."
            
        predictions.append(prediction)
        
        print(f"\nSequence: {seq}")
        print(f" -> 'a' residues: {a_residues}")
        print(f" -> 'd' residues: {d_residues}")
        print(f" -> Analysis: {explanation}")

    print("\n" + "="*30)
    print("Final Predicted Oligomeric States:")
    # Final output needs to be printed clearly.
    result_str = ", ".join(map(str, predictions))
    print(result_str)
    print("="*30)

# The coiled-coil protein sequences provided by the user
protein_sequences = [
    "EIAQALKEIAKALKEIAWALKEIAQALK",
    "EIAALKQEIAALKKENAALKQEIAALKQ",
    "EIAAIKQEIAAIKKEIAAIKWEIAAIKQ",
    "EIQKQLKEIQKQLKEIQWQLKEIQKQLK",
    "EIAQTLKEIAKTLKEIAWTLKEIAQTLK",
]

predict_oligomeric_states(protein_sequences)
<<<D>>>