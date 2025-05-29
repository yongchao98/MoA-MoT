from itertools import permutations

# Possible pieces for each round based on the results
rounds_possibilities = [
    ['D', 'E', 'F', 'G', 'H'],  # Round 1
    ['B', 'C', 'D', 'E', 'F', 'G', 'H'],  # Round 2
    ['A', 'B', 'C', 'D', 'E'],  # Round 3
    ['F', 'G', 'H'],  # Round 4
    ['E', 'F', 'G', 'H']  # Round 5
]

# Check all permutations of opponent's pieces
for perm in permutations('ABCDEFGH', 5):
    if all(perm[i] in rounds_possibilities[i] for i in range(5)):
        print(perm)
        break