# Define the board
board = [
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', 'X', '.', '.', '.'],
    ['X', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', 'O', 'X', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', 'X', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', 'O', 'X', 'X', '.', '.', 'X', '.', '.'],
    ['.', '.', '.', '.', '.', 'X', 'O', 'O', 'X', 'O', '.', '.', '.'],
    ['.', '.', '.', '.', '.', 'X', 'O', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', 'X', '.', '.', '.', '.', '.', '.'],
    ['O', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', 'X', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', 'O', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.']
]

# Function to count liberties of a group
def count_liberties(board, x, y, color):
    visited = set()
    stack = [(x, y)]
    liberties = 0
    while stack:
        cx, cy = stack.pop()
        if (cx, cy) in visited:
            continue
        visited.add((cx, cy))
        for nx, ny in [(cx-1, cy), (cx+1, cy), (cx, cy-1), (cx, cy+1)]:
            if 0 <= nx < len(board) and 0 <= ny < len(board[0]):
                if board[nx][ny] == '.':
                    liberties += 1
                elif board[nx][ny] == color and (nx, ny) not in visited:
                    stack.append((nx, ny))
    return liberties

# Evaluate potential moves
def evaluate_move(board, move_x, move_y):
    # Place the black stone
    board[move_x][move_y] = 'X'
    # Count liberties of the white group at F8 (5, 7)
    liberties = count_liberties(board, 5, 7, 'O')
    # Reset the board
    board[move_x][move_y] = '.'
    return liberties

# Evaluate moves G9 (6, 8) and H7 (7, 7)
liberties_g9 = evaluate_move(board, 6, 8)
liberties_h7 = evaluate_move(board, 7, 7)

# Print the results
print(f"Liberties after G9: {liberties_g9}")
print(f"Liberties after H7: {liberties_h7}")