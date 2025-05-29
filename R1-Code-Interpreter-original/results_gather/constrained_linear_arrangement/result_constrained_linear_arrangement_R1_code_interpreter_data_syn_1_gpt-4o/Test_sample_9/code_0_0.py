from itertools import permutations

# Our cards played in each round
our_cards = ['C', 'A', 'D', 'B']

# Possible opponent cards for each round based on draw conditions
possible_opponent_cards = [
    ['C', 'D'],  # Round 1
    ['A', 'B'],  # Round 2
    ['D', 'A'],  # Round 3
    ['B', 'C']   # Round 4
]

# Check all permutations of opponent's cards
for perm in permutations(['A', 'B', 'C', 'D']):
    if all(perm[i] in possible_opponent_cards[i] for i in range(4)):
        print(list(perm))
        break