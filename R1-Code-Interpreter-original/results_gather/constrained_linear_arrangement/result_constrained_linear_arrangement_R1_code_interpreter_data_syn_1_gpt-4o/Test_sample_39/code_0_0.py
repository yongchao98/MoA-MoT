from itertools import permutations

# Possible pieces for each round based on the analysis
round1_options = ['A', 'B', 'C', 'D', 'E', 'F']
round2_options = ['B', 'C', 'D', 'E', 'F', 'G']
round3_options = ['A', 'B', 'C', 'D', 'E']

# Generate all permutations of opponent pieces
all_permutations = permutations(['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'], 3)

# Find a valid permutation that matches the options for each round
for perm in all_permutations:
    if perm[0] in round1_options and perm[1] in round2_options and perm[2] in round3_options:
        print(perm)
        break