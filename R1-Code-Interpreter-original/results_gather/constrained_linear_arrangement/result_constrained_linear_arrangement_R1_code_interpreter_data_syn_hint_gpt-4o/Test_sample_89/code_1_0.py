# Possible pieces for the opponent
pieces = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']

# Constraints based on each round
constraints = [
    ('A', 'A'),  # Round 1: A vs A (draw)
    ('F', ['B', 'C', 'D', 'E']),  # Round 2: F vs B/C/D/E (loss)
    ('D', ['B', 'C']),  # Round 3: D vs B/C (loss)
    ('C', ['D', 'E', 'F', 'G', 'H']),  # Round 4: C vs D/E/F/G/H (win)
    ('H', ['A', 'B', 'C', 'D', 'E', 'F', 'G'])  # Round 5: H vs A/B/C/D/E/F/G (loss)
]

# Function to check if a permutation satisfies the constraints
def satisfies_constraints(permutation):
    used_pieces = set()
    for i, (my_piece, opponent_options) in enumerate(constraints):
        opponent_piece = permutation[i]
        if opponent_piece in used_pieces:
            return False
        if isinstance(opponent_options, list):
            if opponent_piece not in opponent_options:
                return False
        else:
            if opponent_piece != opponent_options:
                return False
        used_pieces.add(opponent_piece)
    return True

# Find a valid permutation
for perm in permutations(pieces):
    if satisfies_constraints(perm):
        print(perm)
        break