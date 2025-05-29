# Define the pieces and their ranks
pieces = {'A': 8, 'B': 7, 'C': 6, 'D': 5, 'E': 4, 'F': 3, 'G': 2, 'H': 1}

# Our pieces and results
our_pieces = ['F', 'B', 'A', 'H']
results = ['loss', 'win', 'loss', 'loss']

# Opponent's pieces (to be deduced)
opponent_pieces = ['B', 'C', 'H', 'D']

# Verify each round
for i in range(4):
    our_piece = our_pieces[i]
    opponent_piece = opponent_pieces[i]
    result = results[i]
    
    if result == 'win':
        assert pieces[our_piece] > pieces[opponent_piece], f"Round {i+1} failed"
    elif result == 'loss':
        assert pieces[our_piece] < pieces[opponent_piece], f"Round {i+1} failed"
    else:
        assert pieces[our_piece] == pieces[opponent_piece], f"Round {i+1} failed"

print(opponent_pieces)