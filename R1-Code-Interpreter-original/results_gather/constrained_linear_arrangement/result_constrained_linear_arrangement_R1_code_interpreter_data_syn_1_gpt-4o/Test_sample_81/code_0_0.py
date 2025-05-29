# Define the rules
generation = {'A': 'B', 'B': 'C', 'C': 'D', 'D': 'E', 'E': 'A'}
overcoming = {'A': 'C', 'C': 'E', 'E': 'B', 'B': 'D', 'D': 'A'}

# Our cards and results
our_cards = ['D', 'E', 'C', 'B', 'A']
results = ['win', 'draw', 'draw', 'win', 'draw']

# Deduce opponent's cards
opponent_cards = []

for i, result in enumerate(results):
    our_card = our_cards[i]
    if result == 'win':
        # We win if our card overcomes the opponent's card
        for card, overcome in overcoming.items():
            if overcome == our_card:
                opponent_cards.append(card)
                break
    elif result == 'draw':
        # Draw means the opponent played the same card
        opponent_cards.append(our_card)

print(opponent_cards)