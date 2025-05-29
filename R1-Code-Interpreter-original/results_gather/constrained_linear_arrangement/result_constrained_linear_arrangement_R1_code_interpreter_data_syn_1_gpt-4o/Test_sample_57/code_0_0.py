from itertools import permutations

# Our pieces and results
our_pieces = ['E', 'D', 'F', 'A']
results = ['loss', 'win', 'loss', 'loss']

# Possible opponent pieces for each round based on the rules
possible_opponent_pieces = [
    ['A', 'B', 'C', 'D'],  # Round 1: E, loss
    ['E', 'F', 'G', 'H'],  # Round 2: D, win
    ['A', 'B', 'C', 'D', 'E'],  # Round 3: F, loss
    ['H']  # Round 4: A, loss
]

# Check all permutations of opponent pieces
for perm in permutations('ABCDEFGH', 4):
    if (perm[0] in possible_opponent_pieces[0] and
        perm[1] in possible_opponent_pieces[1] and
        perm[2] in possible_opponent_pieces[2] and
        perm[3] in possible_opponent_pieces[3]):
        print(perm)
        break