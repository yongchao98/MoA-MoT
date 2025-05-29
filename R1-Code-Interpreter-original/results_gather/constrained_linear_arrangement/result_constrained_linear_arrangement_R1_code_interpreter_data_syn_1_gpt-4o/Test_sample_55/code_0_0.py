from itertools import permutations

# Possible pieces for each round
round1_options = ['G', 'H']
round2_options = ['B', 'C', 'D', 'E', 'F', 'G', 'H']
round3_options = ['C']
round4_options = ['B', 'C', 'D', 'E', 'F', 'G']

# Generate all permutations of opponent's pieces
all_pieces = ['B', 'C', 'D', 'E', 'F', 'G', 'H']
valid_combinations = []

for perm in permutations(all_pieces, 4):
    if (perm[0] in round1_options and
        perm[1] in round2_options and
        perm[2] in round3_options and
        perm[3] in round4_options):
        valid_combinations.append(perm)

# Print the valid combinations
print(valid_combinations)