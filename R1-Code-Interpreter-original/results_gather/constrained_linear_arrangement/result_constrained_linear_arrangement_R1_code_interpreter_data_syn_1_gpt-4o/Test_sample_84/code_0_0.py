from itertools import permutations

# Our pieces and results
our_pieces = ['D', 'B', 'E', 'G', 'A']
results = ['win', 'win', 'win', 'loss', 'win']

# Possible opponent pieces for each round based on our results
possible_opponent_pieces = [
    ['E', 'F', 'G', 'H'],  # Round 1: D, win
    ['C', 'D', 'E', 'F', 'G', 'H'],  # Round 2: B, win
    ['F', 'G', 'H'],  # Round 3: E, win
    ['A', 'B', 'C', 'D', 'E'],  # Round 4: G, loss
    ['B', 'C', 'D', 'E', 'F', 'G', 'H']  # Round 5: A, win
]

# All opponent pieces
all_opponent_pieces = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']

# Find a valid permutation of opponent pieces
for perm in permutations(all_opponent_pieces, 5):
    if all(perm[i] in possible_opponent_pieces[i] for i in range(5)):
        print(list(perm))
        break