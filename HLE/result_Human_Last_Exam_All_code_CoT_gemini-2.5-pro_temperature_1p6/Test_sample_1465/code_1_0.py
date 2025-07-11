import sys

# Define a function to predict the oligomeric state based on key residues
def predict_coiled_coil_state(sequence):
    """
    Predicts the oligomeric state of a coiled-coil sequence based on
    simplified rules focusing on residues at the 'a' and 'd' positions.
    Assumes a 28-residue sequence with g-a-b-c-d-e-f phasing.
    """
    # 'a' positions are at index 1, 8, 15, 22
    # 'd' positions are at index 4, 11, 18, 25
    a_positions = [1, 8, 15, 22]
    d_positions = [4, 11, 18, 25]

    a_residues = [sequence[i] for i in a_positions]
    d_residues = [sequence[i] for i in d_positions]

    # Rule 1: Asn (N) at an 'a' position strongly favors a dimer.
    if 'N' in a_residues:
        return 2  # Dimer

    # Rule 2: Thr (T) at a 'd' position strongly favors a tetramer.
    if 'T' in d_residues:
        return 4  # Tetramer

    # Rule 3: Otherwise, default to trimer, a common stable state.
    return 3 # Trimer

# The list of protein sequences to analyze
sequences = [
    "EIAQALKEIAKALKEIAWALKEIAQALK",
    "EIAALKQEIAALKKENAALKQEIAALKQ",
    "EIAAIKQEIAAIKKEIAAIKWEIAAIKQ",
    "EIQKQLKEIQKQLKEIQWQLKEIQKQLK",
    "EIAQTLKEIAKTLKEIAWTLKEIAQTLK"
]

# A list to store the results
results = []

# Analyze each sequence and store the result
for seq in sequences:
    state = predict_coiled_coil_state(seq)
    results.append(state)

# Print the final result in the requested format
# Use sys.stdout.write to prevent print's default newline and space for this specific format
output_string = ", ".join(map(str, results))
print(f"The predicted oligomeric states are: {output_string}.")
