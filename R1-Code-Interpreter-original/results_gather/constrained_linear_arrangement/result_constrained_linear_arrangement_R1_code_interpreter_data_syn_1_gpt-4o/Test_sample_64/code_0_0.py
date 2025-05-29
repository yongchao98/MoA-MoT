# Define the hierarchy of pieces
hierarchy = {'A': 8, 'B': 7, 'C': 6, 'D': 5, 'E': 4, 'F': 3, 'G': 2, 'H': 1}

# Your pieces and results
your_pieces = ['H', 'F', 'E', 'C']
results = ['loss', 'loss', 'loss', 'loss']

# Deduce opponent's pieces
opponent_pieces = []

# Available opponent pieces
available_pieces = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']

for i, your_piece in enumerate(your_pieces):
    # Find the opponent's piece that is larger than your piece and not used yet
    for opp_piece in available_pieces:
        if hierarchy[opp_piece] > hierarchy[your_piece]:
            opponent_pieces.append(opp_piece)
            available_pieces.remove(opp_piece)
            break

print(opponent_pieces)