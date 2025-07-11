def predict_coiled_coil_state(sequence):
    """
    Predicts the oligomeric state of a coiled-coil sequence based on established rules.

    The prediction relies on analyzing the amino acids at specific positions
    of the heptad repeat (abcdefg).

    Args:
        sequence (str): The amino acid sequence of the protein.

    Returns:
        int: The predicted oligomeric state (e.g., 2 for dimer, 3 for trimer).
    """
    # Extract residues at key positions of the heptad repeat.
    # We assume the pattern is consistent throughout the sequence.
    a_positions = [sequence[i] for i in range(1, len(sequence), 7)]
    d_positions = [sequence[i] for i in range(4, len(sequence), 7)]
    e_positions = [sequence[i] for i in range(5, len(sequence), 7)]
    c_positions = [sequence[i] for i in range(3, len(sequence), 7)]

    # Rule 1: Trimer (3-mer) is strongly indicated by an Asparagine (N) in the 'a' position.
    if 'N' in a_positions:
        return 3

    # Rule 2: Tetramer (4-mer) is strongly indicated by Isoleucine (I) at both 'a' and 'd' positions.
    if all(p == 'I' for p in a_positions) and all(p == 'I' for p in d_positions):
        return 4

    # Rule 3: For the canonical Isoleucine ('a') / Leucine ('d') core, we check other positions
    # to distinguish between dimers and higher-order oligomers.
    if all(p == 'I' for p in a_positions) and all(p == 'L' for p in d_positions):
        # Rule 3a: A Threonine (T) at the 'e' position can promote a rare heptameric (7-mer) state.
        if 'T' in e_positions:
            return 7
        # Rule 3b: Glutamine (Q) at both 'c' and 'e' positions can promote a pentameric (5-mer) state.
        if 'Q' in e_positions and 'Q' in c_positions:
            return 5
        # Rule 3c: Otherwise, this core defaults to a dimer (2-mer).
        return 2

    # Fallback if no rules match
    return "Unknown"

# The list of protein sequences to analyze
sequences = [
    "EIAQALKEIAKALKEIAWALKEIAQALK",
    "EIAALKQEIAALKKENAALKQEIAALKQ",
    "EIAAIKQEIAAIKKEIAAIKWEIAAIKQ",
    "EIQKQLKEIQKQLKEIQWQLKEIQKQLK",
    "EIAQTLKEIAKTLKEIAWTLKEIAQTLK",
]

# Predict the state for each sequence and store the results
results = []
for seq in sequences:
    state = predict_coiled_coil_state(seq)
    results.append(state)

# Print the final results in the requested format
print("The predicted oligomeric states for the sequences are:")
# The prompt asks to "output each number in the final equation!".
# We will interpret this as printing the final list of numbers.
print(f"{results[0]}, {results[1]}, {results[2]}, {results[3]}, {results[4]}.")