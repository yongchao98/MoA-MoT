from itertools import permutations

# Possible pieces for each round based on the results
round_possibilities = [
    ['C', 'D', 'E', 'F', 'G', 'H'],  # Round 1
    ['A', 'B', 'C', 'D', 'E', 'F'],  # Round 2
    ['B', 'C', 'D', 'E', 'F', 'G'],  # Round 3
    ['A', 'B', 'C', 'D'],            # Round 4
    ['D', 'E', 'F', 'G', 'H']        # Round 5
]

# All possible pieces
all_pieces = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']

# Find a permutation of opponent's pieces that matches the constraints
for perm in permutations(all_pieces, 5):
    if all(perm[i] in round_possibilities[i] for i in range(5)):
        print(perm)
        break