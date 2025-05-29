# Define the initial positions of the knights
white_knights = {'C1': 'w', 'D2': 'w'}
black_knights = {'A2': 'B', 'D3': 'B'}

# Define the target positions for the knights
target_white_knights = {'A2': 'w', 'D3': 'w'}
target_black_knights = {'C1': 'B', 'D2': 'B'}

# Define the moves
moves = [
    ("w", "C1", "B3"),
    ("B", "A2", "C1"),
    ("w", "D2", "A2"),
    ("B", "D3", "B2"),
    ("w", "B3", "D3"),
    ("B", "B2", "D2")
]

# Function to check if the final positions match the target
def check_final_positions(moves, white_knights, black_knights, target_white_knights, target_black_knights):
    # Apply the moves
    for move in moves:
        color, from_pos, to_pos = move
        if color == 'w':
            white_knights[to_pos] = white_knights.pop(from_pos)
        else:
            black_knights[to_pos] = black_knights.pop(from_pos)
    
    # Check if the final positions match the target
    return white_knights == target_white_knights and black_knights == target_black_knights

# Check if the moves achieve the objective
result = check_final_positions(moves, white_knights, black_knights, target_white_knights, target_black_knights)
print("Yes" if result else "No")