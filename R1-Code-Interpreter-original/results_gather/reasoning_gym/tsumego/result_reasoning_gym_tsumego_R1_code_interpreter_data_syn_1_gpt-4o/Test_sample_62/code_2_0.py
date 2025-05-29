def count_liberties(board, x, y):
    """Count the liberties of a stone at position (x, y)."""
    visited = set()
    stack = [(x, y)]
    liberties = 0
    color = board[y][x]

    while stack:
        cx, cy = stack.pop()
        if (cx, cy) in visited:
            continue
        visited.add((cx, cy))

        for nx, ny in [(cx-1, cy), (cx+1, cy), (cx, cy-1), (cx, cy+1)]:
            if 0 <= nx < len(board[0]) and 0 <= ny < len(board):
                if board[ny][nx] == '.':
                    liberties += 1
                elif board[ny][nx] == color and (nx, ny) not in visited:
                    stack.append((nx, ny))

    return liberties

def simulate_move(board, move, color):
    """Simulate a move and return the new board state."""
    x, y = move
    new_board = [list(row) for row in board]
    new_board[y][x] = color
    return new_board

def find_best_move(board):
    """Find the best move for Black to capture White stones."""
    moves = [(6, 8), (7, 8), (12, 8)]  # G8, H8, M8
    best_move = None
    max_captured = 0

    for move in moves:
        new_board = simulate_move(board, move, 'X')
        captured = 0

        # Check the group at G9, H9
        if count_liberties(new_board, 7, 8) == 0:
            captured += 2

        # Check the group at F6, G6, H6
        if count_liberties(new_board, 5, 5) == 0:
            captured += 3

        # Check the isolated stone at M8
        if count_liberties(new_board, 12, 7) == 0:
            captured += 1

        if captured > max_captured:
            max_captured = captured
            best_move = move

    return best_move

# Initial board setup
board = [
    "...........X.",
    ".............",
    ".O...........",
    ".............",
    ".OX....OOX...",
    "...........O.",
    "...........XX",
    ".....O....OOO",
    "..........XO.",
    "..........X..",
    ".............",
    ".............",
    "............."
]

best_move = find_best_move(board)
if best_move:
    columns = "ABCDEFGHIJKLM"
    result = f"<<<{columns[best_move[0]]}{13 - best_move[1]}>>>"
else:
    result = "No valid move found"
print(result)