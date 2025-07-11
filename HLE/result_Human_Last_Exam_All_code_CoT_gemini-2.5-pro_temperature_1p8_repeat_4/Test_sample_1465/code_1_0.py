def predict_oligomeric_state(sequence):
    """
    Predicts the oligomeric state of a coiled-coil sequence based on
    the residues at the core 'a' and 'd' positions.
    For these designed sequences, the phasing is g-a-b-c-d-e-f, and the
    primary determinant is the residue at the first 'd' position (index 4).
    """
    # The residue at the 'a' position is consistently Isoleucine (I).
    # The key determining residue is at the first 'd' position.
    d_residue = sequence[4]

    # Rules are based on empirical data from model coiled-coil studies.
    if d_residue == 'A':
        # (a=I, d=A) pair is a known trimer motif.
        return 3
    elif d_residue == 'L':
        # (a=I, d=L) pair is the canonical dimer signature.
        return 2
    elif d_residue == 'I':
        # (a=I, d=I) pair strongly favors trimers.
        return 3
    elif d_residue == 'Q':
        # A polar (Q) residue in the core 'd' position drives trimerization.
        return 3
    elif d_residue == 'T':
        # A polar (T) residue at 'd' is known to form a tetrameric assembly.
        return 4
    else:
        return "Unknown"

# List of coiled-coil sequences to analyze
sequences = [
    "EIAQALKEIAKALKEIAWALKEIAQALK",
    "EIAALKQEIAALKKENAALKQEIAALKQ",
    "EIAAIKQEIAAIKKEIAAIKWEIAAIKQ",
    "EIQKQLKEIQKQLKEIQWQLKEIQKQLK",
    "EIAQTLKEIAKTLKEIAWTLKEIAQTLK",
]

# Store the final results
final_predictions = []

print("Analyzing coiled-coil sequences:")
for i, seq in enumerate(sequences):
    state = predict_oligomeric_state(seq)
    final_predictions.append(state)
    print(f"Sequence {i+1} ({seq[:7]}...): The predicted oligomeric state is {state}")

# The final "equation" is the sequence of predicted states.
print("\nThe final equation of predicted states is:")
# Output each number as requested
equation_str = ", ".join(map(str, final_predictions))
print(equation_str)