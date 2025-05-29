# Define the cards and their properties
cards = ['A', 'B', 'C', 'D', 'E']
generation = {'A': 'B', 'B': 'C', 'C': 'D', 'D': 'E', 'E': 'A'}
overcoming = {'A': 'C', 'C': 'E', 'E': 'B', 'B': 'D', 'D': 'A'}

# Our cards and results
our_cards = ['B', 'E', 'D', 'C']
results = ['draw', 'win', 'draw', 'draw']

# Deduce opponent's cards
opponent_cards = []

for i, result in enumerate(results):
    if result == 'draw':
        opponent_cards.append(our_cards[i])
    elif result == 'win':
        # Find the card that our card overcomes
        opponent_cards.append(overcoming[our_cards[i]])

print(opponent_cards)