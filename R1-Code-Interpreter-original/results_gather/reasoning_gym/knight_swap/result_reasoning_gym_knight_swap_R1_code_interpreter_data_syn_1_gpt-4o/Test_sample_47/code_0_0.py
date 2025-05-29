# Define the initial positions of the knights
white_knights = {'A2': 'w', 'D1': 'w'}
black_knights = {'A1': 'B', 'B2': 'B'}

# Define the target positions
target_white = {'A1', 'B2'}
target_black = {'A2', 'D1'}

# Define the possible moves for each knight
knight_moves = {
    'A1': ['B3', 'C2'],
    'A2': ['C3', 'C1'],
    'B2': ['D3', 'D1'],
    'D1': ['C3', 'B2']
}

# Function to check if the current configuration matches the target
def is_target_reached(white_knights, black_knights):
    return set(white_knights.keys()) == target_white and set(black_knights.keys()) == target_black

# Simulate the moves
moves = []
current_turn = 'w'  # White moves first

while not is_target_reached(white_knights, black_knights):
    if current_turn == 'w':
        # Try to move a white knight
        for pos, knight in white_knights.items():
            for move in knight_moves[pos]:
                if move not in white_knights and move not in black_knights:
                    # Make the move
                    white_knights[move] = white_knights.pop(pos)
                    moves.append(f"w,{pos},{move}")
                    current_turn = 'B'
                    break
            if current_turn == 'B':
                break
    else:
        # Try to move a black knight
        for pos, knight in black_knights.items():
            for move in knight_moves[pos]:
                if move not in white_knights and move not in black_knights:
                    # Make the move
                    black_knights[move] = black_knights.pop(pos)
                    moves.append(f"B,{pos},{move}")
                    current_turn = 'w'
                    break
            if current_turn == 'w':
                break

# Output the result
if is_target_reached(white_knights, black_knights):
    print(moves)
else:
    print("No")