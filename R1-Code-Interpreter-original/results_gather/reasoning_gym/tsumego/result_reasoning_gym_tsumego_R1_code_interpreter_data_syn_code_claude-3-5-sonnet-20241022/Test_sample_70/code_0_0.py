def print_board(board):
    print("   A B C D E F G H I J")
    for i in range(10):
        row = 10 - i
        print(f"{row:2d}", end=" ")
        for j in range(10):
            print(board[i][j], end=" ")
        print()

def is_valid_move(x, y):
    return 0 <= x < 10 and 0 <= y < 10

def get_liberties(board, x, y, color, visited):
    if not is_valid_move(x, y) or (x, y) in visited:
        return 0
    
    if board[x][y] == '.':
        return 1
    
    if board[x][y] != color:
        return 0
    
    visited.add((x, y))
    liberties = 0
    for dx, dy in [(0, 1), (1, 0), (0, -1), (-1, 0)]:
        liberties += get_liberties(board, x + dx, y + dy, color, visited)
    return liberties

# Initialize the board
board = [
    ['.', '.', '.', '.', '.', '.', '.', '.', 'O', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', 'O', '.', '.', '.', '.', '.', '.', '.', 'O'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', 'X'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', 'O', 'O', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', 'X', 'X', '.'],
    ['.', '.', '.', '.', '.', '.', 'X', 'O', 'O', 'X'],
    ['.', '.', 'O', 'O', '.', '.', 'X', 'O', '.', '.'],
    ['.', '.', 'X', '.', '.', 'X', 'O', 'X', '.', '.']
]

# Check key points
key_points = [(6, 6), (6, 7)]  # G4, H4 coordinates
for x, y in key_points:
    liberties_before = get_liberties(board, x, y, 'O', set())
    print(f"Point ({chr(65+y)}{10-x}): White group has {liberties_before} liberties")