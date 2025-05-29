from itertools import permutations

# Possible pieces for each round based on the analysis
round_possibilities = [
    ['D', 'F', 'G', 'H'],  # Round 1
    ['B', 'C', 'D', 'E', 'F'],  # Round 2
    ['E'],  # Round 3
    ['B', 'C'],  # Round 4
    ['C', 'D', 'F', 'G', 'H']  # Round 5
]

# Generate all permutations of opponent's pieces
all_pieces = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
used_pieces = ['C', 'H', 'E', 'D', 'B']  # Pieces used by the player
remaining_pieces = [p for p in all_pieces if p not in used_pieces]

# Check all permutations of the remaining pieces
for perm in permutations(remaining_pieces, 5):
    if (perm[0] in round_possibilities[0] and
        perm[1] in round_possibilities[1] and
        perm[2] in round_possibilities[2] and
        perm[3] in round_possibilities[3] and
        perm[4] in round_possibilities[4]):
        print(list(perm))
        break