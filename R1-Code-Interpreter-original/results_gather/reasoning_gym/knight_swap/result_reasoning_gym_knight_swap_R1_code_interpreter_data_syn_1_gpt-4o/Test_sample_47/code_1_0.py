# Define the initial positions of the knights
white_knights = {'A2': 'w', 'D1': 'w'}
black_knights = {'A1': 'B', 'B2': 'B'}

# Define the target positions
target_white = {'A1', 'B2'}
target_black = {'A2', 'D1'}

# Function to get possible moves for a knight at a given position
def get_knight_moves(pos):
    col, row = pos
    col_idx = ord(col) - ord('A')
    row_idx = int(row) - 1
    potential_moves = [
        (col_idx + 2, row_idx + 1), (col_idx + 2, row_idx - 1),
        (col_idx - 2, row_idx + 1), (col_idx - 2, row_idx - 1),
        (col_idx + 1, row_idx + 2), (col_idx + 1, row_idx - 2),
        (col_idx - 1, row_idx + 2), (col_idx - 1, row_idx - 2)
    ]
    valid_moves = []
    for c, r in potential_moves:
        if 0 <= c < 4 and 0 <= r < 3:
            valid_moves.append(chr(c + ord('A')) + str(r + 1))
    return valid_moves

# Function to check if the current configuration matches the target
def is_target_reached(white_knights, black_knights):
    return set(white_knights.keys()) == target_white and set(black_knights.keys()) == target_black

# Simulate the moves
moves = []
current_turn = 'w'  # White moves first

while not is_target_reached(white_knights, black_knights):
    if current_turn == 'w':
        # Try to move a white knight
        moved = False
        for pos in list(white_knights.keys()):
            for move in get_knight_moves(pos):
                if move not in white_knights and move not in black_knights:
                    # Make the move
                    white_knights[move] = white_knights.pop(pos)
                    moves.append(f"w,{pos},{move}")
                    current_turn = 'B'
                    moved = True
                    break
            if moved:
                break
    else:
        # Try to move a black knight
        moved = False
        for pos in list(black_knights.keys()):
            for move in get_knight_moves(pos):
                if move not in white_knights and move not in black_knights:
                    # Make the move
                    black_knights[move] = black_knights.pop(pos)
                    moves.append(f"B,{pos},{move}")
                    current_turn = 'w'
                    moved = True
                    break
            if moved:
                break

# Output the result
if is_target_reached(white_knights, black_knights):
    print(moves)
else:
    print("No")