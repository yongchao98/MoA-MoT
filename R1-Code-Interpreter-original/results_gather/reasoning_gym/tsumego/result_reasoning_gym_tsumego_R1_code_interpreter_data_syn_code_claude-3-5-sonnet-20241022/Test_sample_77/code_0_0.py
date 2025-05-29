def print_board(board):
    print("   A B C D E F G H I J K")
    for i in range(11, 0, -1):
        print(f"{i:2d}", end=" ")
        for j in range(11):
            print(board[i-1][j], end=" ")
        print()

def count_liberties(board, x, y, color, visited=None):
    if visited is None:
        visited = set()
    
    if x < 0 or x >= 11 or y < 0 or y >= 11:
        return 0
    
    if board[y][x] != color:
        if board[y][x] == '.':
            return 1
        return 0
    
    if (x, y) in visited:
        return 0
        
    visited.add((x, y))
    liberties = 0
    
    # Check all adjacent positions
    for dx, dy in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
        new_x, new_y = x + dx, y + dy
        liberties += count_liberties(board, new_x, new_y, color, visited)
    
    return liberties

# Initialize the board
board = [
    ['.', 'X', '.', '.', '.', '.', '.', '.', '.', '.', '.'],  # 11
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],  # 10
    ['.', '.', '.', '.', '.', '.', 'O', '.', '.', '.', '.'],  # 9
    ['.', '.', '.', '.', '.', '.', '.', 'X', '.', '.', '.'],  # 8
    ['.', '.', '.', 'X', '.', '.', '.', '.', 'O', '.', '.'],  # 7
    ['.', 'O', 'X', 'O', 'X', '.', '.', '.', '.', '.', '.'],  # 6
    ['.', 'X', 'O', 'O', '.', 'O', '.', '.', '.', '.', '.'],  # 5
    ['.', '.', 'X', 'O', 'O', '.', '.', '.', 'X', '.', '.'],  # 4
    ['.', 'O', '.', 'X', '.', '.', '.', '.', '.', '.', '.'],  # 3
    ['.', '.', '.', 'X', '.', '.', '.', '.', '.', '.', '.'],  # 2
    ['O', '.', 'O', '.', '.', '.', '.', '.', '.', '.', '.']   # 1
]

# Check liberties of key white groups
white_groups = [(1, 5), (2, 4), (3, 4), (4, 3)]  # B6, C5, D5, E4
for x, y in white_groups:
    liberties = count_liberties(board, x, y, 'O')
    print(f"White group at {chr(65+x)}{y+1} has {liberties} liberties")