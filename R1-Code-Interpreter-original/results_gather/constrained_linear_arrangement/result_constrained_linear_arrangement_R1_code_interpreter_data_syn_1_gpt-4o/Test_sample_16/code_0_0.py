# Define the cards and results
your_cards = ['A', 'D', 'C', 'B']
results = ['draw', 'draw', 'win', 'loss']

# Possible opponent cards
opponent_cards = ['A', 'B', 'C', 'D']

# Deduce opponent's cards
opponent_played = []

# Round 1: You played A, result is draw
if results[0] == 'draw':
    opponent_played.append('A' if 'A' in opponent_cards else 'D')

# Round 2: You played D, result is draw
if results[1] == 'draw':
    opponent_played.append('C' if 'C' in opponent_cards else 'D')

# Round 3: You played C, result is win
if results[2] == 'win':
    opponent_played.append('D')

# Round 4: You played B, result is loss
if results[3] == 'loss':
    opponent_played.append('A')

print(opponent_played)