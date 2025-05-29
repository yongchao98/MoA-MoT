from itertools import permutations

# Possible pieces for each round based on the results
round1_possibilities = ['D', 'E', 'F', 'G', 'H']
round2_possibilities = ['B', 'C', 'D', 'E', 'F', 'G']
round3_possibilities = ['C', 'D', 'E', 'F', 'G', 'H']

# Generate all permutations of the opponent's pieces
all_pieces = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
opponent_permutations = permutations(all_pieces, 3)

# Find the valid permutation that matches the conditions
for perm in opponent_permutations:
    if (perm[0] in round1_possibilities and
        perm[1] in round2_possibilities and
        perm[2] in round3_possibilities):
        print(perm)
        break