# Define the pieces and their hierarchy
pieces = {
    'A': 8,  # Elephant
    'B': 7,  # Lion
    'C': 6,  # Tiger
    'D': 5,  # Leopard
    'E': 4,  # Wolf
    'F': 3,  # Dog
    'G': 2,  # Cat
    'H': 1   # Mouse
}

# Results of each round
rounds = [
    ('A', 'loss'),  # Round 1
    ('D', 'win'),   # Round 2
    ('G', 'loss')   # Round 3
]

# Deduce opponent's pieces
opponent_pieces = []

# Round 1: You played 'A' and lost, opponent must have played 'H'
opponent_pieces.append('H')

# Round 2: You played 'D' and won, opponent must have played a piece smaller than 'D'
# Possible pieces: 'E', 'F', 'G', 'H' (but 'H' is already used)
opponent_pieces.append('F')  # Assume 'F' for this round

# Round 3: You played 'G' and lost, opponent must have played a piece larger than 'G'
# Possible pieces: 'A', 'B', 'C', 'D', 'E' (but 'D' is already used by you)
opponent_pieces.append('C')  # Assume 'C' for this round

print(opponent_pieces)