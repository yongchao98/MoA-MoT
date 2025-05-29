from itertools import permutations

# Possible pieces for each round based on the analysis
round_possibilities = [
    ['D', 'C', 'B'],  # Round 1 possibilities
    ['C', 'D', 'E', 'F', 'G', 'H'],  # Round 2 possibilities
    ['B', 'C', 'D', 'E', 'F', 'G'],  # Round 3 possibilities
    ['F', 'G', 'H']  # Round 4 possibilities
]

# All possible opponent pieces
all_pieces = ['B', 'C', 'D', 'E', 'F', 'G', 'H']

# Find a valid permutation of opponent's pieces
for perm in permutations(all_pieces, 4):
    if (perm[0] in round_possibilities[0] and
        perm[1] in round_possibilities[1] and
        perm[2] in round_possibilities[2] and
        perm[3] in round_possibilities[3]):
        print(perm)
        break