from itertools import permutations

# Possible pieces for each round based on the analysis
round1_possibilities = ['B', 'C', 'D', 'E', 'G']
round2_possibilities = ['A', 'B', 'C', 'D', 'E']
round3_piece = 'F'

# Generate all permutations of the opponent's pieces
all_pieces = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
opponent_pieces = [p for p in all_pieces if p != round3_piece]

# Check each permutation to find a valid set of pieces
for perm in permutations(opponent_pieces, 2):
    if perm[0] in round1_possibilities and perm[1] in round2_possibilities:
        print(list(perm) + [round3_piece])
        break