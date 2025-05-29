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
            return {(x,y)}
        return set()
    
    if (x,y) in visited:
        return set()
    
    visited.add((x,y))
    liberties = set()
    
    # Check all adjacent points
    for dx, dy in [(0,1), (0,-1), (1,0), (-1,0)]:
        new_x, new_y = x + dx, y + dy
        liberties.update(get_liberties(board, new_x, new_y, color, visited))
    
    return liberties

# Initialize the board
board = [
    ['.', '.', 'O', '.', '.', '.', '.', 'X', '.'],
    ['.', '.', 'X', '.', '.', '.', 'X', 'O', 'X'],
    ['X', 'X', '.', '.', '.', 'X', 'O', 'O', '.'],
    ['.', '.', '.', 'X', 'X', '.', 'X', 'O', 'X'],
    ['.', '.', '.', '.', '.', '.', '.', 'X', '.'],
    ['.', '.', 'O', '.', '.', 'O', '.', '.', '.'],
    ['X', '.', '.', '.', '.', 'O', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', 'X', '.', '.', '.', '.', '.']
]

# Check liberties of white stones in the critical area
white_group = get_liberties(board, 7, 2, 'O')  # H7
print("Liberties of white group at H7:", white_group)

# Check potential vital points
vital_points = [(5,4), (6,4), (6,5)]  # F5, G5, G6
print("\nVital points to consider:")
for x, y in vital_points:
    print(f"{chr(65+x)}{9-y}: Adjacent to white group")