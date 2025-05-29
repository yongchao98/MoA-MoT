# Define the cards and their properties
cards = ['A', 'B', 'C', 'D', 'E']
generation = {'A': 'B', 'B': 'C', 'C': 'D', 'D': 'E', 'E': 'A'}
overcoming = {'A': 'C', 'C': 'E', 'E': 'B', 'B': 'D', 'D': 'A'}

# Your cards and results
your_cards = ['E', 'D', 'C', 'A']
results = ['draw', 'loss', 'draw', 'draw']

# Determine opponent's cards
opponent_cards = []

for i, your_card in enumerate(your_cards):
    if results[i] == 'draw':
        opponent_cards.append(your_card)
    elif results[i] == 'loss':
        opponent_cards.append(generation[your_card])
    elif results[i] == 'win':
        opponent_cards.append(overcoming[your_card])

print(opponent_cards)