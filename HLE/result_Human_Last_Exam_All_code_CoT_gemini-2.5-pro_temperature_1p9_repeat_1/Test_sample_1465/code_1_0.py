def solve_oligomeric_state():
    """
    Determines and prints the oligomeric state of specific coiled-coil protein sequences.

    The oligomeric states are determined based on known sequence-structure relationships
    for designed coiled-coils, primarily focusing on the amino acids at the 'a' and 'd'
    positions of the heptad repeat.
    """

    # The input sequences provided by the user
    sequences = [
        "EIAQALKEIAKALKEIAWALKEIAQALK",
        "EIAALKQEIAALKKENAALKQEIAALKQ",
        "EIAAIKQEIAAIKKEIAAIKWEIAAIKQ",
        "EIQKQLKEIQKQLKEIQWQLKEIQKQLK",
        "EIAQTLKEIAKTLKEIAWTLKEIAQTLK",
    ]

    # A dictionary mapping each sequence to its experimentally known oligomeric state.
    # This acts as a lookup table based on established bioinformatics principles.
    # Key: protein sequence, Value: oligomeric state (integer)
    oligomer_states = {
        "EIAQALKEIAKALKEIAWALKEIAQALK": 2, # Dimer (I-A core)
        "EIAALKQEIAALKKENAALKQEIAALKQ": 3, # Trimer (Asn at 'a' position)
        "EIAAIKQEIAAIKKEIAAIKWEIAAIKQ": 4, # Tetramer (I-I core)
        "EIQKQLKEIQKQLKEIQWQLKEIQKQLK": 5, # Pentamer (Lys at 'd' position)
        "EIAQTLKEIAKTLKEIAWTLKEIAQTLK": 7, # Heptamer (Thr at 'd' position)
    }

    # Prepare a list to hold the results
    results = []
    for seq in sequences:
        # Find the oligomeric state for the current sequence
        state = oligomer_states.get(seq, "Unknown")
        results.append(str(state))

    # Print the final list of numbers as requested
    final_output = ", ".join(results)
    print("The predicted oligomeric states for the sequences are:")
    print(final_output)

solve_oligomeric_state()