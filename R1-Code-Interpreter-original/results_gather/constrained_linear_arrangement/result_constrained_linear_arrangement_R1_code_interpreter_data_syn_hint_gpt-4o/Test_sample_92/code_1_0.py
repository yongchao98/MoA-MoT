from itertools import permutations

# Possible cards for each round based on the results
round_possibilities = [
    {'B', 'C', 'D', 'E', 'F', 'G'},  # Round 1
    {'B', 'C', 'D', 'E', 'F'},       # Round 2
    {'B', 'C', 'D'},                 # Round 3
    {'C', 'D', 'E', 'F', 'G', 'H'},  # Round 4
    {'E', 'F', 'G', 'H'}             # Round 5
]

# Generate all permutations of the opponent's cards
all_cards = {'B', 'C', 'D', 'E', 'F', 'G', 'H'}
for perm in permutations(all_cards, 5):
    if all(perm[i] in round_possibilities[i] for i in range(5)):
        opponent_cards = list(perm)
        break

print(opponent_cards)