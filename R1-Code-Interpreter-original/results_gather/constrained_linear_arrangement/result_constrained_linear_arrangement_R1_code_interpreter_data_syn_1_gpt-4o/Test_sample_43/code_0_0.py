# Define the cards and results
our_cards = ['A', 'C', 'E', 'D']
results = ['draw', 'draw', 'draw', 'draw']

# Initialize opponent's cards
opponent_cards = []

# Determine opponent's cards based on the results
for i, result in enumerate(results):
    if result == 'draw':
        opponent_cards.append(our_cards[i])

# Print the opponent's cards
print(opponent_cards)