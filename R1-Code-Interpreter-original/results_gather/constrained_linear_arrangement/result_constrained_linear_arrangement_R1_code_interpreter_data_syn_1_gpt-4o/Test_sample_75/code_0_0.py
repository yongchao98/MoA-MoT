# Define the cards and rules
cards = ['A', 'B', 'C', 'D', 'E']
generation = {'A': 'B', 'B': 'C', 'C': 'D', 'D': 'E', 'E': 'A'}
overcoming = {'A': 'C', 'C': 'E', 'E': 'B', 'B': 'D', 'D': 'A'}

# Our cards and results
our_cards = ['B', 'A', 'D', 'E', 'C']
results = ['draw', 'draw', 'draw', 'draw', 'loss']

# Determine opponent's cards
opponent_cards = []

for i, our_card in enumerate(our_cards):
    if results[i] == 'draw':
        opponent_cards.append(our_card)
    elif results[i] == 'loss':
        for card in cards:
            if generation[card] == our_card:
                opponent_cards.append(card)
                break

print(opponent_cards)