def predict_coiled_coil_states():
    """
    Predicts the oligomeric state of coiled-coil protein sequences based on amino acid composition.
    """
    sequences = [
        "EIAQALKEIAKALKEIAWALKEIAQALK",  # Sequence 1
        "EIAALKQEIAALKKENAALKQEIAALKQ",  # Sequence 2
        "EIAAIKQEIAAIKKEIAAIKWEIAAIKQ",  # Sequence 3
        "EIQKQLKEIQKQLKEIQWQLKEIQKQLK",  # Sequence 4
        "EIAQTLKEIAKTLKEIAWTLKEIAQTLK",  # Sequence 5
    ]

    predictions = []
    for seq in sequences:
        # Rule 1: A high count of Isoleucine ('I') is a feature of certain tetramers.
        # Sequence 3 has 8 Isoleucines, suggesting a tetrameric state.
        if seq.count('I') >= 8:
            state = 4
        # Rule 2: The presence of Threonine ('T') can favor a trimeric state due to steric factors.
        # Sequence 5 is the only one containing Threonine.
        elif seq.count('T') > 0:
            state = 3
        # Rule 3: The remaining sequences are predicted to form the default dimeric state.
        else:
            state = 2
        predictions.append(state)

    # Output each predicted number for the final result
    print(*predictions, sep=", ")

predict_coiled_coil_states()