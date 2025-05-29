def print_board(board):
    print("   A B C D E F G H I")
    for i in range(9):
        row = f" {9-i}"
        for j in range(9):
            row += f" {board[i][j]}"
        print(row)

def is_valid_move(x, y):
    return 0 <= x < 9 and 0 <= y < 9

def get_liberties(board, x, y, color, visited=None):
    if visited is None:
        visited = set()
    
    if not is_valid_move(x, y) or (x, y) in visited:
        return set()
    
    if board[x][y] == '.':
        return {(x, y)}
    
    if board[x][y] != color:
        return set()
    
    visited.add((x, y))
    liberties = set()
    
    for dx, dy in [(0, 1), (1, 0), (0, -1), (-1, 0)]:
        new_x, new_y = x + dx, y + dy
        liberties.update(get_liberties(board, new_x, new_y, color, visited))
    
    return liberties

# Initialize the board
board = [
    ['.', '.', '.', '.', '.', '.', '.', '.', 'X'],
    ['.', 'X', 'O', '.', '.', '.', 'O', '.', '.'],
    ['.', '.', '.', '.', 'X', '.', 'X', 'O', 'X'],
    ['.', '.', 'O', 'X', 'O', 'O', 'O', 'X', '.'],
    ['.', '.', 'X', '.', 'X', 'O', 'X', '.', '.'],
    ['.', 'O', 'X', 'X', '.', 'X', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.']
]

# Find white groups and their liberties
for i in range(9):
    for j in range(9):
        if board[i][j] == 'O':
            liberties = get_liberties(board, i, j, 'O')
            if len(liberties) <= 2:  # Groups with 2 or fewer liberties are in danger
                print(f"White group at {chr(65+j)}{9-i} has {len(liberties)} liberties at: {', '.join([chr(65+y) + str(9-x) for x, y in liberties])}")