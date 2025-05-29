# Define the cards and their properties
cards = ['A', 'B', 'C', 'D', 'E']
generation = {'A': 'B', 'B': 'C', 'C': 'D', 'D': 'E', 'E': 'A'}
overcoming = {'A': 'C', 'C': 'E', 'E': 'B', 'B': 'D', 'D': 'A'}

# Your cards and results
your_cards = ['E', 'D', 'C', 'A']
results = ['draw', 'loss', 'draw', 'draw']

# Determine opponent's cards
opponent_cards = [''] * 4

# Round 1: Draw with E
opponent_cards[0] = 'E'

# Round 2: Loss with D, opponent must have played C
opponent_cards[1] = 'C'

# Round 3: Draw with C, opponent must have played B
opponent_cards[2] = 'B'

# Round 4: Draw with A
opponent_cards[3] = 'A'

print(opponent_cards)