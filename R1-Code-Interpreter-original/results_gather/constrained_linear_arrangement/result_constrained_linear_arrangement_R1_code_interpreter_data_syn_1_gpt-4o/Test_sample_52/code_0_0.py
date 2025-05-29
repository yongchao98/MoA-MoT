# Define the rules
generation = {'A': 'B', 'B': 'C', 'C': 'D', 'D': 'E', 'E': 'A'}
overcoming = {'A': 'C', 'C': 'E', 'E': 'B', 'B': 'D', 'D': 'A'}

# Your cards and results
your_cards = ['A', 'B', 'C', 'E']
results = ['win', 'draw', 'win', 'win']

# Determine opponent's cards
opponent_cards = []

for i, your_card in enumerate(your_cards):
    result = results[i]
    if result == 'win':
        # You won, so your card overcomes opponent's card
        for card, overcome in overcoming.items():
            if your_card == card:
                opponent_cards.append(overcome)
                break
    elif result == 'draw':
        # Draw, so opponent played the same card
        opponent_cards.append(your_card)

print(opponent_cards)