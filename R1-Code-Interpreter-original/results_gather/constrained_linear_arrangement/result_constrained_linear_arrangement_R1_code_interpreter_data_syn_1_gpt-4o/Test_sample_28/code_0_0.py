from itertools import permutations

# Possible cards for each round
rounds = [
    ['B', 'D'],  # Round 1
    ['A', 'C'],  # Round 2
    ['D', 'B'],  # Round 3
    ['C', 'A']   # Round 4
]

# Generate all permutations of the cards
all_permutations = permutations(['A', 'B', 'C', 'D'])

# Check each permutation to see if it matches the possible cards for each round
for perm in all_permutations:
    if all(perm[i] in rounds[i] for i in range(4)):
        print(list(perm))
        break