# Define the initial positions of the knights
white_knights = {'A3': 'w', 'C1': 'w'}
black_knights = {'D2': 'B', 'D1': 'B'}

# Define the target positions for the swap
target_white_positions = {'D2', 'D1'}
target_black_positions = {'A3', 'C1'}

# Define the new sequence of moves
moves = [
    ("w", "A3", "B1"),
    ("B", "D2", "C3"),
    ("w", "C1", "A2"),
    ("B", "D1", "B2"),
    ("w", "B1", "D2"),
    ("B", "C3", "A3"),
    ("w", "A2", "D1"),
    ("B", "B2", "C1")
]

# Function to apply moves and check if the swap is successful
def apply_moves(moves, white_knights, black_knights):
    for move in moves:
        color, from_pos, to_pos = move
        if color == 'w':
            if from_pos in white_knights:
                del white_knights[from_pos]
                white_knights[to_pos] = 'w'
        elif color == 'B':
            if from_pos in black_knights:
                del black_knights[from_pos]
                black_knights[to_pos] = 'B'
    
    # Check if the swap is successful
    if set(white_knights.keys()) == target_white_positions and set(black_knights.keys()) == target_black_positions:
        return True
    return False

# Apply the moves and check the result
result = apply_moves(moves, white_knights.copy(), black_knights.copy())
print("Yes" if result else "No")
print(moves if result else [])