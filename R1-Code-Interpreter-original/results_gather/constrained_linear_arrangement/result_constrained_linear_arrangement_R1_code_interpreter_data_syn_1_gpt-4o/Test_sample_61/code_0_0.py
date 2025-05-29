from itertools import permutations

# Possible pieces for each round based on the results
round1_possibilities = {'A', 'B', 'C', 'D', 'E'}
round2_possibilities = {'E', 'F', 'G', 'H'}
round3_possibilities = {'B', 'C', 'D', 'E', 'F', 'G'}
round4_possibilities = {'B', 'C', 'D', 'E', 'F', 'G'}

# Generate all permutations of opponent's pieces
all_pieces = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
possible_permutations = permutations(all_pieces, 4)

# Find the valid permutation that satisfies all conditions
for perm in possible_permutations:
    if (perm[0] in round1_possibilities and
        perm[1] in round2_possibilities and
        perm[2] in round3_possibilities and
        perm[3] in round4_possibilities):
        print(perm)
        break