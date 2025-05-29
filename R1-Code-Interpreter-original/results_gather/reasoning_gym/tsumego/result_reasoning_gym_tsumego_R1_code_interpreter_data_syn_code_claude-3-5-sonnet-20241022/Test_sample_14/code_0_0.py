def print_board(board):
    print("   A B C D E F G H I")
    for i in range(9):
        row = f" {9-i}"
        for j in range(9):
            row += f" {board[i][j]}"
        print(row)

def get_liberties(board, x, y, color, visited=None):
    if visited is None:
        visited = set()
    
    if x < 0 or x >= 9 or y < 0 or y >= 9:
        return set()
    
    if board[y][x] != color:
        if board[y][x] == '.':
            return {(x, y)}
        return set()
    
    if (x, y) in visited:
        return set()
    
    visited.add((x, y))
    liberties = set()
    
    # Check all adjacent points
    for dx, dy in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
        new_x, new_y = x + dx, y + dy
        liberties.update(get_liberties(board, new_x, new_y, color, visited))
        if new_x >= 0 and new_x < 9 and new_y >= 0 and new_y < 9 and board[new_y][new_x] == '.':
            liberties.add((new_x, new_y))
    
    return liberties

# Initialize the board
board = [
    ['.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', 'O', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', 'O'],
    ['.', '.', '.', '.', '.', '.', '.', '.', 'X'],
    ['.', '.', '.', 'X', '.', '.', '.', '.', '.'],
    ['O', '.', 'X', 'O', 'X', '.', '.', '.', 'X'],
    ['O', 'X', 'O', 'O', '.', '.', '.', '.', '.'],
    ['.', 'X', 'X', 'O', 'X', 'X', '.', '.', '.'],
    ['X', 'X', '.', 'X', '.', 'O', '.', '.', '.']
]

# Find white group and its liberties
white_group_liberties = get_liberties(board, 0, 6, 'O')
print("White group liberties:", [(chr(x+65), 9-y) for x, y in white_group_liberties])