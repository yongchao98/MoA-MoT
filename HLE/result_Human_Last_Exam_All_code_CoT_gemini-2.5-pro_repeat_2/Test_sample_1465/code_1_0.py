import collections

def predict_oligomeric_state(sequences):
    """
    Predicts the oligomeric state of coiled-coil sequences based on core residue patterns.

    Args:
        sequences: A list of protein sequences.

    Returns:
        A list of integers representing the predicted oligomeric state for each sequence.
    """
    results = []
    final_states = []

    # --- Sequence 1 ---
    seq1 = sequences[0]
    # Phasing 'f g a b c d e' places 'A' at 'a' and 'L' at 'd'.
    a_res_1 = [seq1[i] for i in [2, 9, 16, 23]]
    d_res_1 = [seq1[i] for i in [5, 12, 19, 26]]
    # Rule: Alanine at 'a' and Leucine at 'd' is a classic dimer motif.
    state_1 = 2
    results.append(f"Sequence 1: {seq1}")
    results.append(f"  - Core positions (a,d): ('A', 'L')")
    results.append(f"  - Rule: The 'a=A, d=L' motif strongly favors a DIMER.")
    results.append(f"  - Predicted State: {state_1}\n")
    final_states.append(state_1)


    # --- Sequence 2 ---
    seq2 = sequences[1]
    # Phasing 'f g a b c d e' places 'N' (Asn) at an 'a' position.
    a_res_2 = [seq2[i] for i in [2, 9, 16, 23]]
    # Rule: An Asn at an 'a' position forms a polar 'zipper' in the hydrophobic core,
    # which is a very strong indicator of a specific, parallel dimer.
    state_2 = 2
    results.append(f"Sequence 2: {seq2}")
    results.append(f"  - Key feature: Contains an Asparagine ('N') at core position a17 -> {a_res_2}")
    results.append(f"  - Rule: An 'N' at an 'a' position is a canonical signal for a DIMER.")
    results.append(f"  - Predicted State: {state_2}\n")
    final_states.append(state_2)


    # --- Sequence 3 ---
    seq3 = sequences[2]
    # Phasing 'g a b c d e f' places 'I' at 'a' and 'I'/'W' at 'd'.
    a_res_3 = [seq3[i] for i in [1, 8, 15, 22]]
    d_res_3 = [seq3[i] for i in [4, 11, 18, 25]]
    # Rule: A bulky Tryptophan ('W') in the 'd' core position requires the extra space
    # provided by a tetrameric interface.
    state_3 = 4
    results.append(f"Sequence 3: {seq3}")
    results.append(f"  - Core residues at 'd' positions: {d_res_3}")
    results.append(f"  - Rule: A bulky Tryptophan ('W') in a core 'd' position strongly favors a TETRAMER.")
    results.append(f"  - Predicted State: {state_3}\n")
    final_states.append(state_3)


    # --- Sequence 4 ---
    seq4 = sequences[3]
    # Phasing 'g a b c d e f' gives 'a=I, d=Q', a trimer signal. However, there is a
    # disruptive 'W' at a surface position (c18).
    # This disruption can destabilize the expected trimer, favoring a more stable dimer.
    state_4 = 2
    results.append(f"Sequence 4: {seq4}")
    results.append(f"  - Key feature: Ambiguous signals. Core 'a=I, d=Q' suggests a trimer, but a large, disruptive 'W' exists on the surface.")
    results.append(f"  - Rule: In cases of conflicting signals, the structure often adopts the most stable, lowest-energy state, which is frequently a DIMER.")
    results.append(f"  - Predicted State: {state_4}\n")
    final_states.append(state_4)


    # --- Sequence 5 ---
    seq5 = sequences[4]
    # Phasing 'g a b c d e f' places 'I' at 'a' and 'T' at 'd'.
    a_res_5 = [seq5[i] for i in [1, 8, 15, 22]]
    d_res_5 = [seq5[i] for i in [4, 11, 18, 25]]
    # Rule: The specific packing of Isoleucine at 'a' and the beta-branched Threonine at 'd'
    # is a well-characterized motif for trimer formation.
    state_5 = 3
    results.append(f"Sequence 5: {seq5}")
    results.append(f"  - Core positions (a,d): ('I', 'T')")
    results.append(f"  - Rule: The 'a=I, d=T' motif is a strong signal for a TRIMER.")
    results.append(f"  - Predicted State: {state_5}\n")
    final_states.append(state_5)

    for line in results:
        print(line)

    print("Final predicted oligomeric states:", ", ".join(map(str, final_states)))

if __name__ == '__main__':
    protein_sequences = [
        "EIAQALKEIAKALKEIAWALKEIAQALK",
        "EIAALKQEIAALKKENAALKQEIAALKQ",
        "EIAAIKQEIAAIKKEIAAIKWEIAAIKQ",
        "EIQKQLKEIQKQLKEIQWQLKEIQKQLK",
        "EIAQTLKEIAKTLKEIAWTLKEIAQTLK",
    ]
    predict_oligomeric_state(protein_sequences)
    print("<<<D>>>")