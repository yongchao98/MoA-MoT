# Define the cards and results
my_cards = ['B', 'A', 'C', 'D']
results = ['draw', 'loss', 'draw', 'loss']

# Define the rules
rules = {
    'A': 'D',  # Stick beats Worm
    'B': 'A',  # Tiger beats Stick
    'C': 'B',  # Chicken beats Tiger
    'D': 'C'   # Worm beats Chicken
}

# Determine opponent's cards
opponent_cards = []

for i, my_card in enumerate(my_cards):
    if results[i] == 'draw':
        opponent_cards.append(my_card)
    elif results[i] == 'loss':
        for card, beats in rules.items():
            if beats == my_card:
                opponent_cards.append(card)
                break

print(opponent_cards)