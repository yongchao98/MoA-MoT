from itertools import permutations

# Possible pieces for the opponent
pieces = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']

# Conditions for each round
conditions = [
    lambda x: x in ['C', 'D', 'E', 'F', 'G', 'H'],  # Round 1: smaller than B
    lambda x: x in ['A', 'B', 'C', 'D', 'E', 'F'],  # Round 2: larger than G
    lambda x: x in ['B', 'C', 'D', 'E', 'F', 'G'],  # Round 3: larger than H
    lambda x: x in ['A', 'B', 'C', 'D', 'E'],       # Round 4: larger than F
    lambda x: x in ['D', 'E', 'F', 'G', 'H']        # Round 5: smaller than C
]

# Check all permutations of 5 pieces
for perm in permutations(pieces, 5):
    if all(cond(perm[i]) for i, cond in enumerate(conditions)):
        print(perm)
        break