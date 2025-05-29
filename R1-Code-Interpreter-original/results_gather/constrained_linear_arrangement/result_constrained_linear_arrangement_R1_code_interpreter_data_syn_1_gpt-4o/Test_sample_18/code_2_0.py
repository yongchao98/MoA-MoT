# Define the cards and results
my_cards = ['B', 'A', 'C', 'D']
results = ['draw', 'loss', 'draw', 'loss']

# Define the rules
rules = {
    'A': 'D',  # Stick is beaten by Worm
    'B': 'B',  # Tiger is beaten by Tiger (draw)
    'C': 'C',  # Chicken is beaten by Chicken (draw)
    'D': 'C'   # Worm is beaten by Chicken
}

# Determine opponent's cards
opponent_cards = []

for i, my_card in enumerate(my_cards):
    if results[i] == 'draw':
        opponent_cards.append(my_card)
    elif results[i] == 'loss':
        opponent_cards.append(rules[my_card])

print(opponent_cards)