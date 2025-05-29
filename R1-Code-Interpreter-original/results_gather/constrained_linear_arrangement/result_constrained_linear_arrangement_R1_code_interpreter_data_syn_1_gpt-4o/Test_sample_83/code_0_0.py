# Define the cards and their properties
cards = {
    'A': {'generates': 'C', 'overcomes': 'B'},
    'B': {'generates': 'D', 'overcomes': 'E'},
    'C': {'generates': 'E', 'overcomes': 'D'},
    'D': {'generates': 'A', 'overcomes': 'C'},
    'E': {'generates': 'B', 'overcomes': 'A'}
}

# Our plays and results
our_plays = ['B', 'E', 'A', 'D', 'C']
results = ['draw', 'draw', 'loss', 'draw', 'draw']

# Deduce opponent's plays
opponent_plays = []

for i, our_card in enumerate(our_plays):
    if results[i] == 'draw':
        opponent_plays.append(our_card)
    elif results[i] == 'loss':
        opponent_plays.append(cards[our_card]['generates'])

print(opponent_plays)