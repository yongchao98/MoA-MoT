# Define the cards and results
my_cards = ['B', 'A', 'C', 'D']
results = ['draw', 'loss', 'draw', 'loss']

# Define the rules
rules = {
    'A': 'B',  # Stick loses to Tiger
    'B': 'C',  # Tiger loses to Chicken
    'C': 'D',  # Chicken loses to Worm
    'D': 'A'   # Worm loses to Stick
}

# Determine opponent's cards
opponent_cards = []

for i, my_card in enumerate(my_cards):
    if results[i] == 'draw':
        opponent_cards.append(my_card)
    elif results[i] == 'loss':
        opponent_cards.append(rules[my_card])

print(opponent_cards)