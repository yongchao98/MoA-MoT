# Define the cards and results
your_cards = ['A', 'B', 'C', 'D', 'E']
results = ['draw', 'draw', 'draw', 'draw', 'draw']

# Determine opponent's cards based on the results
opponent_cards = []

for i, result in enumerate(results):
    if result == 'draw':
        opponent_cards.append(your_cards[i])

print(opponent_cards)