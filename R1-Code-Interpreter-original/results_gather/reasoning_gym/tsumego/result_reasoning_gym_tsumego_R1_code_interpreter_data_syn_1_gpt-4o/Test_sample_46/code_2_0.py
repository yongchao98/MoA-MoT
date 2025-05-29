# Directly check the specific move that seems promising
def check_specific_move(board, move):
    x, y = move
    new_board = [list(row) for row in board]
    new_board[y][x] = 'X'
    white_groups = parse_board(new_board)
    captures = 0
    for group, liberties in white_groups:
        if not liberties:
            captures += len(group)
    return captures

board = [
    "..........",
    "..X.O.....",
    "..O...O...",
    "..XO...X..",
    "..O.X.O...",
    ".XOOOXO...",
    "..XOX.....",
    "...X......",
    "O.........",
    ".........."
]

# Check the move D9
move = (3, 1)  # D9 in 0-indexed coordinates
captures = check_specific_move(board, move)
columns = "ABCDEFGHIJ"
if captures > 0:
    x, y = move
    print(f"<<<{columns[x]}{10-y}>>>")
else:
    print("No capturing move found.")