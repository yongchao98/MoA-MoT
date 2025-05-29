def print_board(board):
    print("   A B C D E F G H I J K")
    for i in range(11, 0, -1):
        row = f"{i:2d} "
        for j in range(11):
            row += board[i-1][j] + " "
        print(row)

def get_liberties(board, x, y, color, visited=None):
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
        liberties.update(get_liberties(board, new_x, new_y, color, visited))
        if new_x >= 0 and new_x < 11 and new_y >= 0 and new_y < 11 and board[new_y][new_x] == '.':
            liberties.add((new_x, new_y))
    
    return liberties

# Initialize the board
board = [
    ['.'] * 11 for _ in range(11)
]

# Set up the position
position = """
. . . . X . X . . . .
. X . X . . . X X . .
X O O O X . . X . . .
X X O X . . . . . . .
. . X . . . . . . O .
. . . . O . . . . . .
. . O . . . . O . . .
. . . . . . . . . . .
. . . X . . . . . . .
. . . . . . . . . . .
. . . . . . . . . . .
"""

# Parse the position into the board
rows = position.strip().split('\n')
for i, row in enumerate(rows):
    cells = row.split()
    for j, cell in enumerate(cells):
        board[i][j] = cell

# Find liberties of the white group at C9
white_liberties = get_liberties(board, 2, 8, 'O')  # C9 coordinates

print("White group liberties:")
for x, y in white_liberties:
    col = chr(65 + x)
    row = 11 - y
    print(f"{col}{row}")