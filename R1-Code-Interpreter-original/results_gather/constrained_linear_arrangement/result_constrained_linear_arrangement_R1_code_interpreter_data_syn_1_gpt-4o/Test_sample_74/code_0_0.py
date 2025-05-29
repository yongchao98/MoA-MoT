# Our cards played in each round
our_cards = ['D', 'E', 'A', 'B', 'C']

# Results of each round
results = ['draw', 'draw', 'draw', 'draw', 'draw']

# Opponent's cards based on the results
opponent_cards = []

# Since all results are draw, opponent's cards must match ours
for i in range(len(our_cards)):
    if results[i] == 'draw':
        opponent_cards.append(our_cards[i])

print(opponent_cards)