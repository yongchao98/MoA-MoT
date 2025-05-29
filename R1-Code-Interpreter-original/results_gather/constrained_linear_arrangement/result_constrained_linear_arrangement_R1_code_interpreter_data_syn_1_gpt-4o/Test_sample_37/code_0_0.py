from itertools import permutations

# Your cards and results
your_cards = ['A', 'C', 'D', 'B']
results = ['loss', 'draw', 'draw', 'win']

# Possible outcomes based on the rules
outcomes = {
    'A': {'win': 'B', 'loss': 'D', 'draw': 'A'},
    'B': {'win': 'C', 'loss': 'A', 'draw': 'B'},
    'C': {'win': 'D', 'loss': 'B', 'draw': 'C'},
    'D': {'win': 'A', 'loss': 'C', 'draw': 'D'}
}

# Check all permutations of opponent's cards
for opponent_cards in permutations(['A', 'B', 'C', 'D']):
    valid = True
    for i in range(4):
        if outcomes[your_cards[i]][results[i]] != opponent_cards[i]:
            valid = False
            break
    if valid:
        print(list(opponent_cards))
        break