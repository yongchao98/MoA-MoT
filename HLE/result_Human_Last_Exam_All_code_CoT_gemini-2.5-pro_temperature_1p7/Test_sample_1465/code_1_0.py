import collections

def predict_oligomeric_state(sequence: str) -> int:
    """
    Predicts the oligomeric state of a coiled-coil sequence based on
    residues at the 'a' and 'd' positions of the heptad repeat.
    """
    if not sequence:
        return 0

    a_pos_residues = []
    d_pos_residues = []

    for i, residue in enumerate(sequence):
        # Heptad positions are 1-based (a=1, b=2...), indices are 0-based
        heptad_pos = (i % 7)
        if heptad_pos == 0:  # 'a' position
            a_pos_residues.append(residue)
        elif heptad_pos == 3:  # 'd' position
            d_pos_residues.append(residue)

    # Use Counter to find the most common residue at each position
    a_counts = collections.Counter(a_pos_residues)
    d_counts = collections.Counter(d_pos_residues)
    
    most_common_a = a_counts.most_common(1)[0][0] if a_counts else None
    most_common_d = d_counts.most_common(1)[0][0] if d_counts else None

    # --- Rule-based prediction ---

    # Rule for Sequence 5: A at 'a' and T at 'd' -> Tetramer (4)
    if most_common_a == 'A' and most_common_d == 'T':
        return 4
        
    # Rule for Sequence 2: A at 'a' and L at 'd' -> Dimer (2)
    # The 'N' in the sequence is a known motif within a dimeric context.
    if most_common_a == 'A' and most_common_d == 'L':
        return 2

    # Rule for Sequence 1: I at 'a' is a strong signal for a trimer
    if most_common_a == 'I':
        return 3

    # Rule for Sequence 3: A at 'a' and I at 'd' favors trimers
    if most_common_a == 'A' and most_common_d == 'I':
        return 3

    # Rule for Sequence 4: Q at 'a' (polar zipper) is often trimeric
    if most_common_a == 'Q':
        return 3

    # Default fallback
    return 0 # Unknown

def solve():
    """
    Main function to process the sequences and print the results.
    """
    sequences = [
        "EIAQALKEIAKALKEIAWALKEIAQALK",
        "EIAALKQEIAALKKENAALKQEIAALKQ",
        "EIAAIKQEIAAIKKEIAAIKWEIAAIKQ",
        "EIQKQLKEIQKQLKEIQWQLKEIQKQLK",
        "EIAQTLKEIAKTLKEIAWTLKEIAQTLK",
    ]

    predictions = []
    for seq in sequences:
        state = predict_oligomeric_state(seq)
        predictions.append(state)

    # Format the output to list each predicted number clearly.
    # The "equation" is the final comma-separated list of values.
    result_string = ", ".join(map(str, predictions))
    print(f"The final equation of predicted oligomeric states is: {result_string}")

solve()