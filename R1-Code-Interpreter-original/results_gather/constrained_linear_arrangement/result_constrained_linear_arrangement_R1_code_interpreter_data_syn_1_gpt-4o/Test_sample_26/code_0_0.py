from itertools import permutations

# Possible cards for each round based on the draw condition
round_possibilities = [
    ['D', 'C'],  # Round 1
    ['C', 'B'],  # Round 2
    ['A', 'D'],  # Round 3
    ['B', 'A']   # Round 4
]

# Generate all permutations of the opponent's cards
all_permutations = permutations(['A', 'B', 'C', 'D'])

# Check each permutation to see if it matches the round possibilities
for perm in all_permutations:
    if (perm[0] in round_possibilities[0] and
        perm[1] in round_possibilities[1] and
        perm[2] in round_possibilities[2] and
        perm[3] in round_possibilities[3]):
        print(list(perm))
        break