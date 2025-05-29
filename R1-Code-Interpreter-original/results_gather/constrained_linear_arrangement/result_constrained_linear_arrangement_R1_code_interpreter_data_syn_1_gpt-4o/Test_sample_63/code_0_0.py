# Define the pieces and their hierarchy
pieces = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']

# Results of each round
rounds = [
    ('G', 'loss'),  # Round 1
    ('D', 'loss'),  # Round 2
    ('C', 'loss'),  # Round 3
    ('A', 'win')    # Round 4
]

# Deduce opponent's pieces
opponent_pieces = ['C', 'A', 'B', 'H']  # Based on reasoning

print(opponent_pieces)