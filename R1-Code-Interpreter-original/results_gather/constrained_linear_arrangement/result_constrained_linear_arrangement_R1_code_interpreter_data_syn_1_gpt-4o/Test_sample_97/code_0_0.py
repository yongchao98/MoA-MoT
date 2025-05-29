from itertools import permutations

# Possible pieces for the opponent
pieces = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']

# Constraints based on the rounds
constraints = [
    ('A', 'A'),  # Round 1: Draw with A
    ('E', ['B', 'C', 'D']),  # Round 2: Loss with E
    ('D', ['E', 'F', 'G', 'H']),  # Round 3: Win with D
    ('H', ['B', 'C', 'D', 'E', 'F', 'G']),  # Round 4: Loss with H
    ('F', ['G', 'H'])  # Round 5: Win with F
]

# Function to check if a permutation satisfies the constraints
def satisfies_constraints(permutation):
    for i, (my_piece, possible_opponent_pieces) in enumerate(constraints):
        if isinstance(possible_opponent_pieces, list):
            if permutation[i] not in possible_opponent_pieces:
                return False
        else:
            if permutation[i] != possible_opponent_pieces:
                return False
    return True

# Find the valid permutation
for perm in permutations(pieces, 5):
    if satisfies_constraints(perm):
        print(list(perm))
        break