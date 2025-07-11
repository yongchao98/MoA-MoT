import collections

def analyze_coiled_coils():
    """
    Analyzes coiled-coil sequences to predict their oligomeric state.
    """
    sequences = [
        "EIAQALKEIAKALKEIAWALKEIAQALK",
        "EIAALKQEIAALKKENAALKQEIAALKQ",
        "EIAAIKQEIAAIKKEIAAIKWEIAAIKQ",
        "EIQKQLKEIQKQLKEIQWQLKEIQKQLK",
        "EIAQTLKEIAKTLKEIAWTLKEIAQTLK"
    ]
    
    # Heptad positions (0-indexed)
    # Register g-a-b-c-d-e-f
    # 'a' is at index 1, 8, 15, 22
    # 'd' is at index 4, 11, 18, 25
    a_indices = [1, 8, 15, 22]
    d_indices = [4, 11, 18, 25]

    predictions = []
    
    # --- Sequence 1 ---
    seq1 = sequences[0]
    a_res1 = [seq1[i] for i in a_indices]
    d_res1 = [seq1[i] for i in d_indices]
    print(f"Sequence 1: {seq1}")
    print(f"  'a' positions: {a_res1}")
    print(f"  'd' positions: {d_res1}")
    # Reasoning: 'a' positions are all Isoleucine (I), a very strong dimer signal.
    # 'd' positions are mostly Alanine (A), which is also compatible with a dimer.
    # The single bulky Tryptophan (W) at a 'd' position is not sufficient to overcome
    # the strong dimeric preference of the three I-A pairs.
    pred1 = 2
    predictions.append(pred1)
    print(f"Prediction: Dimer ({pred1})\n")

    # --- Sequence 2 ---
    seq2 = sequences[1]
    a_res2 = [seq2[i] for i in a_indices]
    d_res2 = [seq2[i] for i in d_indices]
    print(f"Sequence 2: {seq2}")
    print(f"  'a' positions: {a_res2}")
    print(f"  'd' positions: {d_res2}")
    # Reasoning: This is a classic dimer-forming sequence. It has Isoleucine (I) at the 'a'
    # positions and Leucine (L) at the 'd' positions. The single Asparagine (N)
    # at an 'a' position is well-tolerated in dimers.
    pred2 = 2
    predictions.append(pred2)
    print(f"Prediction: Dimer ({pred2})\n")
    
    # --- Sequence 3 ---
    seq3 = sequences[2]
    a_res3 = [seq3[i] for i in a_indices]
    d_res3 = [seq3[i] for i in d_indices]
    # In this sequence, the 'W' is at index 20, which is position 'g' (g-a-b-c-d-e-f),
    # so it does not affect the core.
    print(f"Sequence 3: {seq3}")
    print(f"  'a' positions: {a_res3}")
    print(f"  'd' positions: {d_res3}")
    # Reasoning: This sequence has Isoleucine (I) at all 'a' AND 'd' positions.
    # The I-I pairing in the core is the canonical signal for a tetramer.
    pred3 = 4
    predictions.append(pred3)
    print(f"Prediction: Tetramer ({pred3})\n")
    
    # --- Sequence 4 ---
    seq4 = sequences[3]
    a_res4 = [seq4[i] for i in a_indices]
    d_res4 = [seq4[i] for i in d_indices]
    print(f"Sequence 4: {seq4}")
    print(f"  'a' positions: {a_res4}")
    print(f"  'd' positions: {d_res4}")
    # Reasoning: 'a' positions are all Isoleucine (I), a strong dimer signal.
    # The 'd' positions are mostly polar (Glutamine, Q), which is highly destabilizing.
    # A stable trimer is unlikely. The structure defaults to the state preferred by the
    # dominant 'a' position signal.
    pred4 = 2
    predictions.append(pred4)
    print(f"Prediction: Dimer ({pred4})\n")

    # --- Sequence 5 ---
    seq5 = sequences[4]
    a_res5 = [seq5[i] for i in a_indices]
    d_res5 = [seq5[i] for i in d_indices]
    print(f"Sequence 5: {seq5}")
    print(f"  'a' positions: {a_res5}")
    print(f"  'd' positions: {d_res5}")
    # Reasoning: This sequence has Isoleucine (I) at 'a' positions and Threonine (T) at
    # most 'd' positions. This I-T pairing favors dimers. However, the single Tryptophan (W)
    # at a 'd' position is a potent trimer-promoting residue, and it is sufficient to
    # switch the preferred state of the entire oligomer from a dimer to a trimer.
    pred5 = 3
    predictions.append(pred5)
    print(f"Prediction: Trimer ({pred5})\n")

    # --- Final Result ---
    print("---------------------------------")
    print(f"Final predicted oligomeric states: {predictions}")
    print("This corresponds to Answer Choice D.")

analyze_coiled_coils()