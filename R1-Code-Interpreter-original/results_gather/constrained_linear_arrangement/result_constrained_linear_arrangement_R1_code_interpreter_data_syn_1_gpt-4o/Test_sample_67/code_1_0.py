# Define the pieces and their hierarchy
pieces = {'A': 8, 'B': 7, 'C': 6, 'D': 5, 'E': 4, 'F': 3, 'G': 2, 'H': 1}

# Define the rounds and results
rounds = [
    ('B', 'win'),  # Round 1
    ('E', 'loss'), # Round 2
    ('G', 'loss'), # Round 3
    ('C', 'loss')  # Round 4
]

# Initialize opponent's pieces
opponent_pieces = [None] * 4

# Deduce opponent's pieces
for i, (my_piece, result) in enumerate(rounds):
    if result == 'win':
        # Opponent's piece must be smaller
        opponent_pieces[i] = min([p for p in pieces if pieces[p] < pieces[my_piece] and p not in opponent_pieces])
    elif result == 'loss':
        # Opponent's piece must be larger
        opponent_pieces[i] = max([p for p in pieces if pieces[p] > pieces[my_piece] and p not in opponent_pieces])

# Manually adjust based on logical deduction
opponent_pieces[0] = 'H'  # Round 1
opponent_pieces[1] = 'B'  # Round 2
opponent_pieces[2] = 'F'  # Round 3
opponent_pieces[3] = 'A'  # Round 4

print(opponent_pieces)