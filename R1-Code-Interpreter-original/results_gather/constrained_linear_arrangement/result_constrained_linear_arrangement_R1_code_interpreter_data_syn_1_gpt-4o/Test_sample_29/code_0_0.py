from itertools import permutations

# Possible cards
cards = ['A', 'B', 'C', 'D']

# Check all permutations of the opponent's cards
for perm in permutations(cards):
    if (perm[0] in ['B', 'A'] and  # Round 1
        perm[1] in ['C', 'D'] and  # Round 2
        perm[2] in ['C', 'D'] and  # Round 3
        perm[3] in ['A', 'B']):    # Round 4
        print(perm)
        break