from itertools import permutations

# Possible pieces for each round
rounds = [
    ['G', 'F', 'E'],  # Round 1
    ['A'],            # Round 2
    ['D', 'C', 'B'],  # Round 3
    ['E', 'F', 'G', 'H'],  # Round 4
    ['D', 'E', 'F', 'G', 'H']  # Round 5
]

# Generate all permutations of opponent's pieces
all_pieces = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
used_pieces = set()

# Check each permutation
for perm in permutations(all_pieces, 5):
    if (perm[0] in rounds[0] and
        perm[1] in rounds[1] and
        perm[2] in rounds[2] and
        perm[3] in rounds[3] and
        perm[4] in rounds[4]):
        used_pieces = perm
        break

print(used_pieces)