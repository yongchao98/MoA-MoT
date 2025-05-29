# Define the hierarchy of pieces
pieces = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']

# Results of each round
results = ['loss', 'win', 'loss', 'win', 'win']

# Our pieces played in each round
our_pieces = ['E', 'G', 'D', 'B', 'A']

# Opponent's possible pieces for each round
opponent_pieces = [[], [], [], [], []]

# Deduce opponent's pieces based on the results
for i, result in enumerate(results):
    if result == 'loss':
        # Opponent's piece must be higher than ours
        opponent_pieces[i] = pieces[:pieces.index(our_pieces[i])]
    elif result == 'win':
        # Opponent's piece must be lower than ours
        opponent_pieces[i] = pieces[pieces.index(our_pieces[i]) + 1:]

# Known piece from Round 2
opponent_pieces[1] = ['H']

# Deduce remaining pieces
used_pieces = set(opponent_pieces[1])
for i in range(len(opponent_pieces)):
    if len(opponent_pieces[i]) > 1:
        for piece in opponent_pieces[i]:
            if piece not in used_pieces:
                opponent_pieces[i] = [piece]
                used_pieces.add(piece)
                break

# Print the deduced opponent's pieces
print([opponent_pieces[i][0] for i in range(len(opponent_pieces))])