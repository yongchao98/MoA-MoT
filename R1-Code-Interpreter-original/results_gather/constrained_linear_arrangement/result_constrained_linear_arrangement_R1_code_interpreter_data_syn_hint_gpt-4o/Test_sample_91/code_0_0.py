from itertools import permutations

# Possible pieces for each round based on the analysis
round1_options = ['A', 'B', 'C']
round3_options = ['B', 'C']
round4_options = ['C', 'E', 'F', 'G']
round5_options = ['B', 'C', 'E', 'F', 'G']

# All possible pieces
all_pieces = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']

# Function to check if a permutation is valid
def is_valid_permutation(perm):
    return (perm[0] in round1_options and
            perm[2] in round3_options and
            perm[3] in round4_options and
            perm[4] in round5_options)

# Find the valid permutation
for perm in permutations(all_pieces, 5):
    if perm[1] == 'H' and is_valid_permutation(perm):
        print(perm)
        break