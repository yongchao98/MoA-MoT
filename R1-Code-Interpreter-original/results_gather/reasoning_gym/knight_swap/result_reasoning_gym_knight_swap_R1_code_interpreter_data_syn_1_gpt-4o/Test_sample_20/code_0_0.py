# Define the initial positions of the knights
white_knights = {'A2': 'w', 'D2': 'w'}
black_knights = {'C3': 'B', 'C2': 'B'}

# Define the target positions for the knights
target_white_knights = {'C3': 'w', 'C2': 'w'}
target_black_knights = {'A2': 'B', 'D2': 'B'}

# Define the moves
moves = [
    ("w", "A2", "C3"),
    ("B", "C3", "A2"),
    ("w", "D2", "C2"),
    ("B", "C2", "D2")
]

# Function to check if the moves lead to the target configuration
def check_moves(moves, white_knights, black_knights, target_white_knights, target_black_knights):
    for move in moves:
        color, from_pos, to_pos = move
        if color == 'w':
            if from_pos in white_knights:
                del white_knights[from_pos]
                white_knights[to_pos] = 'w'
            else:
                return False
        elif color == 'B':
            if from_pos in black_knights:
                del black_knights[from_pos]
                black_knights[to_pos] = 'B'
            else:
                return False

    return white_knights == target_white_knights and black_knights == target_black_knights

# Check if the sequence of moves is valid
is_valid = check_moves(moves, white_knights.copy(), black_knights.copy(), target_white_knights, target_black_knights)
print("Yes" if is_valid else "No")
if is_valid:
    print(moves)