from itertools import permutations

# Possible pieces for the opponent
pieces = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']

# Constraints based on each round
constraints = [
    ['G', 'F', 'E', 'D', 'C', 'B', 'A'],  # Round 1: H, loss
    ['F', 'E', 'D', 'C', 'B', 'A'],       # Round 2: G, loss
    ['D', 'C', 'B', 'A'],                 # Round 3: E, loss
    ['C', 'D', 'E', 'F', 'G', 'H'],       # Round 4: B, win
    ['E', 'F', 'G', 'H']                  # Round 5: D, win
]

# Check all permutations of opponent's pieces
for perm in permutations(pieces, 5):
    if (perm[0] in constraints[0] and
        perm[1] in constraints[1] and
        perm[2] in constraints[2] and
        perm[3] in constraints[3] and
        perm[4] in constraints[4]):
        print(perm)
        break