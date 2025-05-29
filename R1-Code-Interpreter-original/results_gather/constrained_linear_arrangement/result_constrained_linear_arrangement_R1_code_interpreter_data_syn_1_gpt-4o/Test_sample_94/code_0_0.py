from itertools import permutations

# Possible pieces for each round based on the analysis
round_possibilities = [
    ['A', 'B', 'C', 'E', 'F'],  # Round 1
    ['A', 'B', 'C', 'E'],       # Round 2
    ['C', 'E', 'F', 'G'],       # Round 3
    ['D'],                      # Round 4
    ['H']                       # Round 5
]

# All pieces the opponent can use
all_pieces = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']

# Find a permutation of opponent's pieces that matches the constraints
for perm in permutations(all_pieces):
    if (perm[0] in round_possibilities[0] and
        perm[1] in round_possibilities[1] and
        perm[2] in round_possibilities[2] and
        perm[3] == 'D' and
        perm[4] == 'H'):
        opponent_pieces = perm[:5]
        break

print(opponent_pieces)