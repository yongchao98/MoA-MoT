def print_board(board):
    print("   A B C D E F G H I J K")
    for i in range(11, 0, -1):
        row = f"{i:2d} "
        for j in range(11):
            row += board[i-1][j] + " "
        print(row)

def get_group_liberties(board, x, y, color, visited=None):
    if visited is None:
        visited = set()
    
    if x < 0 or x >= 11 or y < 0 or y >= 11:
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
        if 0 <= new_x < 11 and 0 <= new_y < 11:
            if board[new_y][new_x] == '.':
                liberties.add((new_x, new_y))
            elif board[new_y][new_x] == color:
                liberties.update(get_group_liberties(board, new_x, new_y, color, visited))

    return liberties

# Initialize the board
board = [
    ['.', '.', '.', '.', 'X', '.', 'X', '.', '.', '.', '.'],  # 11
    ['.', 'X', '.', 'X', '.', '.', '.', 'X', 'X', '.', '.'],  # 10
    ['X', 'O', 'O', 'O', 'X', '.', '.', 'X', '.', '.', '.'],  # 9
    ['X', 'X', 'O', 'X', '.', '.', '.', '.', '.', '.', '.'],  # 8
    ['.', '.', 'X', '.', '.', '.', '.', '.', '.', 'O', '.'],  # 7
    ['.', '.', '.', '.', 'O', '.', '.', '.', '.', '.', '.'],  # 6
    ['.', '.', 'O', '.', '.', '.', '.', 'O', '.', '.', '.'],  # 5
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],  # 4
    ['.', '.', '.', 'X', '.', '.', '.', '.', '.', '.', '.'],  # 3
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],  # 2
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.']   # 1
]

# Find liberties of the white group starting at C9 (2,8)
white_group_liberties = get_group_liberties(board, 2, 2, 'O')

print("White group liberties:")
for x, y in white_group_liberties:
    col = chr(65 + x)
    row = 11 - y
    print(f"{col}{row}")

# Also print the board for verification
print("\nCurrent board position:")
print_board(board)