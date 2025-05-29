def is_surrounded(board, x, y, color):
    if x < 0 or x >= 9 or y < 0 or y >= 9:
        return True
    if board[y][x] == '.':
        return False
    if board[y][x] != color:
        return True
    return True

def analyze_move(board, move_x, move_y):
    # Create a copy of the board with the move
    new_board = [row[:] for row in board]
    new_board[move_y][move_x] = 'X'
    
    # Check liberties of surrounding white groups
    white_liberties = 0
    for y in range(max(0, move_y-1), min(9, move_y+2)):
        for x in range(max(0, move_x-1), min(9, move_x+2)):
            if new_board[y][x] == 'O':
                if not (is_surrounded(new_board, x+1, y, 'O') and
                       is_surrounded(new_board, x-1, y, 'O') and
                       is_surrounded(new_board, x, y+1, 'O') and
                       is_surrounded(new_board, x, y-1, 'O')):
                    white_liberties += 1
    
    return white_liberties

# Initialize the board
board = [
    ['.', 'X', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['O', '.', '.', 'O', '.', '.', '.', '.', '.'],
    ['X', '.', '.', '.', '.', '.', '.', 'X', '.'],
    ['X', '.', '.', '.', '.', '.', 'X', 'O', '.'],
    ['.', '.', '.', 'X', '.', 'X', '.', '.', '.'],
    ['.', '.', '.', '.', 'O', 'O', 'X', '.', '.'],
    ['.', '.', '.', 'X', 'O', 'O', '.', 'X', '.'],
    ['.', '.', '.', '.', 'X', 'O', 'X', '.', '.']
]

# Check E3 move
liberties = analyze_move(board, 4, 6)  # E3 coordinates
print(f"White liberties after Black plays at E3: {liberties}")