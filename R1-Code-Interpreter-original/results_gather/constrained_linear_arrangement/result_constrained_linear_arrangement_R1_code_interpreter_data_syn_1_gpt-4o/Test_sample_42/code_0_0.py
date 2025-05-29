# Define the rules for generation and overcoming
generation = {'A': 'B', 'B': 'C', 'C': 'D', 'D': 'E', 'E': 'A'}
overcoming = {'A': 'C', 'C': 'E', 'E': 'B', 'B': 'D', 'D': 'A'}

# Our cards and results
our_cards = ['A', 'D', 'E', 'C']
results = ['loss', 'win', 'draw', 'loss']

# Determine opponent's cards
opponent_cards = []

for i, our_card in enumerate(our_cards):
    result = results[i]
    if result == 'loss':
        # We lose if our card generates the opponent's card
        opponent_card = generation[our_card]
    elif result == 'win':
        # We win if our card overcomes the opponent's card
        opponent_card = overcoming[our_card]
    elif result == 'draw':
        # Draw if both cards are the same
        opponent_card = our_card
    opponent_cards.append(opponent_card)

print(opponent_cards)