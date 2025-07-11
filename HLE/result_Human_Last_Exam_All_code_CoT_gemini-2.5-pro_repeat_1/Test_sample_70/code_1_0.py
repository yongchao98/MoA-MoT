import numpy as np

def get_matrix_rank(matrix):
    """Calculates the rank of a given matrix."""
    return np.linalg.matrix_rank(np.array(matrix))

# --- Define Flag Matrices ---
# The matrices represent the geometric structure of the flags.
# We assign unique numbers to distinct color/pattern regions to ensure maximal rank.

# Denmark's flag has a Scandinavian cross structure, which creates two
# independent row/column patterns.
# Let Red=1, White=2. Matrix: [[R,W,R], [W,W,W], [R,W,R]]
denmark_matrix = [[1, 2, 1], [2, 2, 2], [1, 2, 1]]

# A selection of African flags with varied structures.
african_flags = {
    # Rank 1: Simple horizontal/vertical stripes.
    'Nigeria': [[1, 2, 1], [1, 2, 1], [1, 2, 1]],  # Vertical tricolor (Green=1, White=2)
    'Gabon': [[1, 1, 1], [2, 2, 2], [3, 3, 3]],      # Horizontal tricolor (Green=1, Yellow=2, Blue=3)

    # Rank 2: Hoist bar + horizontal stripes.
    'Benin': [[1, 2, 2], [1, 3, 3]],  # Vertical Green(1) bar, horizontal Yellow(2)/Red(3) stripes.
    'Madagascar': [[1, 2, 2], [1, 3, 3]], # Vertical White(1) bar, horizontal Red(2)/Green(3) stripes.
    # Rank 2: Hoist bar with emblem + horizontal stripes.
    'Guinea-Bissau': [[1, 2, 2], [3, 4, 4]], # Red(1,3) bar w/ Star, Yellow(2)/Green(4) stripes.

    # Rank 2: Canton + striped field.
    'Liberia': [[1, 2, 2], [3, 3, 3], [2, 2, 2]], # Blue/White Canton(1), Red(2)/White(3) stripes.
    'Togo': [[1, 2, 2], [3, 3, 3], [2, 2, 2]], # Red/White Canton(1), Green(2)/Yellow(3) stripes.

    # Rank 2: Emblem on a plain field.
    'Somalia': [[1, 1, 1], [1, 2, 1], [1, 1, 1]], # White(2) star on a Blue(1) field.
    'Tunisia': [[1, 1, 1], [1, 2, 1], [1, 1, 1]], # White(2) disk/emblem on a Red(1) field.

    # Rank 3: Complex geometric designs.
    'South Africa': [[1, 2, 3], [4, 4, 4], [5, 6, 4]], # Y-shape design with 6 colors.
    'Eritrea': [[1, 1, 2], [1, 2, 3], [2, 3, 3]], # Hoist triangle dividing two other triangles.
}

# 1. Calculate the rank of the Danish flag.
target_rank = get_matrix_rank(denmark_matrix)
print(f"The rank of the flag of Denmark is: {target_rank}")
print("-" * 30)

# 2. Find African flags with the same rank.
matching_flags = []
for country, matrix in african_flags.items():
    rank = get_matrix_rank(matrix)
    if rank == target_rank:
        matching_flags.append(country)

# 3. Print the result.
print(f"African nations with flags of rank {target_rank} are:")
for country in sorted(matching_flags):
    print(f"- {country}")
