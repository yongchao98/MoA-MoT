from itertools import permutations

# Possible pieces for each round based on the results
round1_possibilities = ['B', 'C', 'D', 'E', 'F', 'G']
round2_possibilities = ['D', 'E', 'F', 'G', 'H']
round3_possibilities = ['A', 'B', 'C']
round4_possibilities = ['A', 'B', 'C', 'D', 'E', 'F']

# Generate all permutations of 4 pieces from the opponent's set
opponent_pieces = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
for perm in permutations(opponent_pieces, 4):
    if (perm[0] in round1_possibilities and
        perm[1] in round2_possibilities and
        perm[2] in round3_possibilities and
        perm[3] in round4_possibilities):
        print(perm)
        break