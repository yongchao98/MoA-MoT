from itertools import permutations

# Possible pieces for each round based on the analysis
rounds_possibilities = [
    ['A', 'B', 'C', 'D'],  # Round 1: F, loss
    ['A', 'B', 'C', 'D'],  # Round 2: E, loss
    ['E', 'F', 'G', 'H'],  # Round 3: D, win
    ['C'],                 # Round 4: C, draw
    ['B', 'D', 'E', 'F', 'G', 'H']  # Round 5: A, win
]

# All possible pieces
all_pieces = set('ABCDEFGH')

# Function to check if a permutation is valid
def is_valid_permutation(perm):
    # Check if the permutation matches the possibilities for each round
    for i, piece in enumerate(perm):
        if piece not in rounds_possibilities[i]:
            return False
    return True

# Generate all permutations of the opponent's pieces
for perm in permutations(all_pieces, 5):
    if is_valid_permutation(perm):
        print(perm)
        break