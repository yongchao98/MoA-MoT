from itertools import permutations

# Possible pieces for each round based on the analysis
round_possibilities = [
    ['A', 'B', 'C', 'D', 'E'],  # Round 1: F, loss
    ['D'],                      # Round 2: D, draw
    ['F', 'G', 'H'],            # Round 3: E, win
    ['G']                       # Round 4: H, loss
]

# All pieces available to the opponent
all_pieces = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']

# Function to check if a permutation of pieces is valid
def is_valid_permutation(permutation):
    for i, piece in enumerate(permutation):
        if piece not in round_possibilities[i]:
            return False
    return True

# Generate all permutations of the opponent's pieces
for perm in permutations(all_pieces, 4):
    if is_valid_permutation(perm):
        print(perm)
        break