# Define the cards and their properties
cards = ['A', 'B', 'C', 'D', 'E']
generation = {'A': 'B', 'B': 'C', 'C': 'D', 'D': 'E', 'E': 'A'}
overcoming = {'A': 'C', 'C': 'E', 'E': 'B', 'B': 'D', 'D': 'A'}

# Our plays and results
our_plays = ['B', 'A', 'E', 'C', 'D']
results = ['draw', 'loss', 'draw', 'win', 'draw']

# Deduce opponent's plays
opponent_plays = []

for i, our_card in enumerate(our_plays):
    result = results[i]
    if result == 'draw':
        opponent_plays.append(our_card)
    elif result == 'loss':
        opponent_plays.append(generation[our_card])
    elif result == 'win':
        for card, overcomes in overcoming.items():
            if overcomes == our_card:
                opponent_plays.append(card)
                break

print(opponent_plays)