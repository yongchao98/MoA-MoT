# Define the cards and results
my_cards = ['C', 'A', 'B', 'D']
results = ['loss', 'draw', 'loss', 'draw']

# Define the rules
rules = {
    'A': 'B',  # Stick beats Tiger
    'B': 'C',  # Tiger beats Chicken
    'C': 'D',  # Chicken beats Worm
    'D': 'A'   # Worm beats Stick
}

# Determine opponent's cards
opponent_cards = []

for i, result in enumerate(results):
    my_card = my_cards[i]
    if result == 'draw':
        opponent_cards.append(my_card)
    elif result == 'loss':
        # Find the card that beats my card
        for card, beats in rules.items():
            if beats == my_card:
                opponent_cards.append(card)
                break

print(opponent_cards)