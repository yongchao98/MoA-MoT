# Define the rounds and results
rounds = [
    ('H', 'loss'),
    ('A', 'win'),
    ('F', 'draw'),
    ('G', 'draw')
]

# Define the hierarchy and special rule
hierarchy = {'A': 8, 'B': 7, 'C': 6, 'D': 5, 'E': 4, 'F': 3, 'G': 2, 'H': 1}
special_rule = {'H': 'A'}  # Mouse can capture Elephant

# Deduce opponent's pieces
opponent_pieces = []

for my_piece, result in rounds:
    if result == 'loss':
        # Find the piece that can beat my piece
        if my_piece == 'H':
            opponent_pieces.append('G')  # Mouse loses to Cat
        else:
            for piece, rank in hierarchy.items():
                if rank > hierarchy[my_piece]:
                    opponent_pieces.append(piece)
                    break
    elif result == 'win':
        # Find the largest piece that my piece can beat
        for piece, rank in hierarchy.items():
            if rank < hierarchy[my_piece]:
                opponent_pieces.append(piece)
                break
    elif result == 'draw':
        # Same piece as mine
        opponent_pieces.append(my_piece)

print(opponent_pieces)