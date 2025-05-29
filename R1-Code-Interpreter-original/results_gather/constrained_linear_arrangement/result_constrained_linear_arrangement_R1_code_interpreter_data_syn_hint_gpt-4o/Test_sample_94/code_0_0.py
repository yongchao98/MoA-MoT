from itertools import permutations

# Possible opponent pieces
opponent_pieces = ['A', 'B', 'C', 'E', 'F', 'G']

# Check all permutations of the opponent's pieces
for perm in permutations(opponent_pieces, 3):
    round1, round2, round3 = perm
    if (round1 in ['A', 'B', 'C', 'E', 'F'] and  # Round 1: larger than G
        round2 in ['A', 'B', 'C', 'E'] and      # Round 2: larger than F
        round3 in ['C', 'E', 'F', 'G']):        # Round 3: smaller than B
        print([round1, round2, 'D', round3, 'H'])