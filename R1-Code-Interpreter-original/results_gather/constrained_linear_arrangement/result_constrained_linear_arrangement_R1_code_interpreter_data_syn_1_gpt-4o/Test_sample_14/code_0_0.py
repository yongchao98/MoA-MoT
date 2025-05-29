from itertools import permutations

# Possible cards for each round based on the results
rounds = [
    ('A', 'D'),  # Round 1: Draw with A, so opponent could have played A or D
    ('B', 'A'),  # Round 2: Draw with B, so opponent could have played B or A
    ('B',),      # Round 3: Loss with C, so opponent must have played B
    ('C',)       # Round 4: Win with D, so opponent must have played C
]

# All possible cards opponent can play
cards = ['A', 'B', 'C', 'D']

# Find a permutation of cards that matches the possible outcomes
for perm in permutations(cards):
    if (perm[0] in rounds[0] and
        perm[1] in rounds[1] and
        perm[2] in rounds[2] and
        perm[3] in rounds[3]):
        print(list(perm))
        break