from itertools import permutations

# Possible pieces for each round based on the results
round1_possibilities = ['A', 'B', 'C']
round2_possibilities = ['C', 'D', 'E', 'F', 'G', 'H']
round3_possibilities = ['D', 'E', 'F', 'G', 'H']
round4_possibilities = ['B', 'C', 'D', 'E', 'F', 'G']

# All pieces
all_pieces = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']

# Find a valid permutation of opponent's pieces
for perm in permutations(all_pieces, 4):
    if (perm[0] in round1_possibilities and
        perm[1] in round2_possibilities and
        perm[2] in round3_possibilities and
        perm[3] in round4_possibilities):
        print(perm)
        break