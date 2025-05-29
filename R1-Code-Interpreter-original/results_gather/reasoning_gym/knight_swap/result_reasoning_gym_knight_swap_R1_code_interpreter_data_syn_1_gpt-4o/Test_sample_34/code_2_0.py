# Define the initial positions of the knights
white_knights = {'B1': 'w', 'D2': 'w'}
black_knights = {'A3': 'B', 'A2': 'B'}
empty_squares = ['C3', 'C2']

# Define the new moves
new_moves = [
    ("w", "B1", "C3"),
    ("B", "A3", "C1"),
    ("w", "D2", "C2"),
    ("B", "A2", "C3"),
    ("w", "C2", "A2"),
    ("B", "C1", "D2"),
    ("w", "C3", "A3"),
    ("B", "C3", "B1")
]

# Function to apply moves and check if the swap is successful
def apply_moves(moves, white_knights, black_knights, empty_squares):
    for move in moves:
        color, from_pos, to_pos = move
        if color == 'w':
            if from_pos in white_knights and to_pos in empty_squares:
                white_knights[to_pos] = white_knights.pop(from_pos)
                empty_squares.remove(to_pos)
                empty_squares.append(from_pos)
            else:
                return "No"
        elif color == 'B':
            if from_pos in black_knights and to_pos in empty_squares:
                black_knights[to_pos] = black_knights.pop(from_pos)
                empty_squares.remove(to_pos)
                empty_squares.append(from_pos)
            else:
                return "No"
    
    # Check if the final positions match the objective
    if white_knights == {'A3': 'w', 'A2': 'w'} and black_knights == {'B1': 'B', 'D2': 'B'}:
        return moves
    else:
        return "No"

# Apply the new moves and print the result
result = apply_moves(new_moves, white_knights, black_knights, empty_squares)
print(result)