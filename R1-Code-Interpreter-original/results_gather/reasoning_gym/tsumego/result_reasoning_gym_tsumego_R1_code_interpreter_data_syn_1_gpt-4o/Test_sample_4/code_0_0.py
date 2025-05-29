# Define the board
board = [
    ['O', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', 'O', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', 'X', '.', '.', '.', '.', '.', '.', 'X', '.', '.'],
    ['.', '.', '.', '.', '.', '.', 'O', '.', 'X', '.', '.'],
    ['.', '.', 'X', '.', '.', '.', '.', 'O', 'O', 'X', '.'],
    ['.', '.', '.', 'O', '.', '.', 'X', 'O', 'O', '.', 'X'],
    ['.', '.', '.', '.', '.', '.', '.', 'X', 'O', 'X', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', 'X', '.', '.'],
    ['.', '.', 'X', '.', '.', '.', 'X', '.', '.', '.', '.']
]

# Check the liberties of the group at D4 and E4
def check_liberties(board, group):
    liberties = set()
    for (r, c) in group:
        for dr, dc in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
            nr, nc = r + dr, c + dc
            if 0 <= nr < len(board) and 0 <= nc < len(board[0]) and board[nr][nc] == '.':
                liberties.add((nr, nc))
    return liberties

# Group at D4 and E4
group = [(7, 3), (7, 4)]
liberties = check_liberties(board, group)

# Determine the move
if len(liberties) == 1:
    move = liberties.pop()
    move_coord = (chr(move[1] + ord('A')), 11 - move[0])
    print(f"<<<{move_coord[0]}{move_coord[1]}>>>")
else:
    print("No capture possible")
