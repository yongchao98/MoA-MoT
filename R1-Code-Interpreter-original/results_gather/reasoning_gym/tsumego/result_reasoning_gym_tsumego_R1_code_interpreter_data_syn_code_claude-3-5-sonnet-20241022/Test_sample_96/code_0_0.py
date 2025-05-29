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
    ['.', '.', 'X', '.', 'X', '.', '.', '.', '.', '.'],  # 10
    ['X', 'X', 'O', 'O', 'O', 'X', '.', '.', '.', '.'],  # 9
    ['.', '.', 'X', 'O', 'X', '.', '.', '.', '.', '.'],  # 8
    ['.', '.', '.', 'X', '.', '.', '.', '.', '.', '.'],  # 7
    ['.', '.', '.', 'O', '.', '.', '.', '.', '.', '.'],  # 6
    ['.', 'O', '.', '.', '.', '.', 'O', 'X', 'X', '.'],  # 5
    ['.', '.', '.', '.', '.', '.', '.', '.', 'X', '.'],  # 4
    ['O', '.', 'O', '.', '.', '.', '.', '.', '.', '.'],  # 3
    ['.', '.', '.', '.', '.', 'O', '.', '.', '.', '.'],  # 2
    ['.', '.', '.', '.', '.', '.', '.', '.', 'X', '.']   # 1
]

# Check liberties of white group
white_group_liberties = set()
for i in range(10):
    for j in range(10):
        if board[i][j] == 'O':
            liberties = get_liberties(board, i, j, 'O', set())
            if liberties < 3:  # Looking for vulnerable groups
                print(f"White group at {chr(65+j)}{10-i} has {liberties} liberties")

print("\nCurrent board position:")
print_board(board)