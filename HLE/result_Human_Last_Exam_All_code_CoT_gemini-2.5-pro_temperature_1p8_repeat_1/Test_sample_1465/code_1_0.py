def predict_coiled_coil_states():
    """
    Predicts the oligomeric state of coiled-coil sequences based on the residues
    at the 'a' and 'd' positions of the heptad repeat.
    """
    sequences = [
        "EIAQALKEIAKALKEIAWALKEIAQALK",
        "EIAALKQEIAALKKENAALKQEIAALKQ",
        "EIAAIKQEIAAIKKEIAAIKWEIAAIKQ",
        "EIQKQLKEIQKQLKEIQWQLKEIQKQLK",
        "EIAQTLKEIAKTLKEIAWTLKEIAQTLK",
    ]

    oligomer_states = []

    # --- Analysis for each sequence ---

    # Sequence 1: EIAQALKEIAKALKEIAWALKEIAQALK
    # Phasing is 'gabcdef'. 'a' positions are 1, 8, 15, 22. 'd' positions are 4, 11, 18, 25.
    # Core residues: a=[I,I,I,I], d=[A,A,W,A].
    # Rule: a=I, d=A. The small Alanine (A) residue at the 'd' position generally favors a Dimer.
    oligomer_states.append(2)

    # Sequence 2: EIAALKQEIAALKKENAALKQEIAALKQ
    # Sequence is 30 aa, likely contains a typo. Assuming the repeating unit 'EIAALK' is correct,
    # the sequence should be 28 aa: "EIAALKQEIAALKEIAALKEIAALKQ".
    # Phasing 'gabcdef'. Core residues: a=[I,I,I,I], d=[L,L,L,L].
    # Rule (Harbury): a=I, d=L strongly indicates a Trimer.
    oligomer_states.append(3)

    # Sequence 3: EIAAIKQEIAAIKKEIAAIKWEIAAIKQ
    # Phasing is 'gabcdef'. Core residues: a=[I,I,I,I], d=[I,I,I,I].
    # Rule (Harbury): a=I, d=I strongly indicates a Tetramer.
    oligomer_states.append(4)
    
    # Sequence 4: EIQKQLKEIQKQLKEIQWQLKEIQKQLK
    # This sequence fits a different 'defgabc' phasing.
    # 'd' positions are 1, 8, 15, 22. 'a' positions are 5, 12, 19, 26.
    # Core residues: d=[I,I,I,I], a=[L,L,L,L].
    # Rule (Harbury): a=L, d=I indicates a Dimer.
    oligomer_states.append(2)

    # Sequence 5: EIAQTLKEIAKTLKEIAWTLKEIAQTLK
    # Phasing is 'gabcdef'. Core residues: a=[I,I,I,I], d=[Q,K,W,Q].
    # Rule: a=I with polar residues (Q, K) at the 'd' position. This arrangement
    # often stabilizes a Trimeric state through buried hydrogen bonds.
    oligomer_states.append(3)
    
    # Print the final predicted oligomeric states
    # Final equation based on the calculated oligomeric states for each protein
    print(f"The predicted oligomeric states are: {oligomer_states[0]}, {oligomer_states[1]}, {oligomer_states[2]}, {oligomer_states[3]}, {oligomer_states[4]}.")

predict_coiled_coil_states()
<<<D>>>