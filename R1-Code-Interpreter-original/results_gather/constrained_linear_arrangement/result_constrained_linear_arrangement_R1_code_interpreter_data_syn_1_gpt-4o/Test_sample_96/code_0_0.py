from itertools import permutations

# Possible pieces for each round based on the analysis
round_possibilities = [
    ['F', 'G', 'H'],  # Round 1
    ['C'],            # Round 2
    ['E', 'F', 'G', 'H'],  # Round 3
    ['D', 'E', 'F', 'G', 'H'],  # Round 4
    ['A', 'B', 'D', 'E']   # Round 5
]

# Generate all permutations of opponent's pieces
all_permutations = permutations('ABCDEFGH', 5)

# Find the valid permutation that matches the round possibilities
for perm in all_permutations:
    if (perm[0] in round_possibilities[0] and
        perm[1] in round_possibilities[1] and
        perm[2] in round_possibilities[2] and
        perm[3] in round_possibilities[3] and
        perm[4] in round_possibilities[4]):
        print(perm)
        break