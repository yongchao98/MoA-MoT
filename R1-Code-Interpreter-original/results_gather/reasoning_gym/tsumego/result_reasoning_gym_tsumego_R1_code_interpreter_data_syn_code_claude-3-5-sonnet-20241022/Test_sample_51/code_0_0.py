def print_board(board):
    print("   A B C D E F G H I J K L")
    for i in range(12, 0, -1):
        print(f"{i:2d}", end=" ")
        for j in range(12):
            print(board[i-1][j], end=" ")
        print()

def count_captures(board, x, y):
    if board[y][x] != '.':
        return 0
    
    # Temporarily place black stone
    board[y][x] = 'X'
    
    # Check for captures
    captures = 0
    directions = [(0,1), (1,0), (0,-1), (-1,0)]
    visited = set()
    
    def has_liberty(x, y, visited):
        if x < 0 or x >= 12 or y < 0 or y >= 12:
            return False
        if board[y][x] == '.':
            return True
        if board[y][x] == 'X':
            return False
        if (x,y) in visited:
            return False
        
        visited.add((x,y))
        for dx, dy in directions:
            if has_liberty(x+dx, y+dy, visited):
                return True
        return False

    # Check each adjacent position for white groups without liberties
    for dx, dy in directions:
        nx, ny = x + dx, y + dy
        if 0 <= nx < 12 and 0 <= ny < 12 and board[ny][nx] == 'O':
            visited = set()
            if not has_liberty(nx, ny, visited):
                captures += len(visited)
    
    # Remove test stone
    board[y][x] = '.'
    return captures

# Initialize board
board = [['.' for _ in range(12)] for _ in range(12)]

# Set up the position
position = """
12 O . . . . . . X . . . .
11 X . . . . . . . . . . .
10 O . . . . . . . . . . .
 9 . . . X . . . O . . . .
 8 . . . . . . . . . . . .
 7 . . . . O . . X . . . .
 6 . . . . . . . . . . . .
 5 . X . . . . . . . . . .
 4 . . . . . . X . . . . .
 3 . . . . . . X X . . . X
 2 . . . . . X O O X . . .
 1 . . . O . X O . . . . .
"""

# Parse position into board
rows = position.strip().split('\n')
for row in rows:
    parts = row.split()
    y = int(parts[0]) - 1
    for x, char in enumerate(parts[1:]):
        board[y][x] = char

# Find best capturing move
max_captures = 0
best_move = None

for y in range(12):
    for x in range(12):
        captures = count_captures(board, x, y)
        if captures > max_captures:
            max_captures = captures
            best_move = (x, y)

if best_move:
    col = chr(best_move[0] + ord('A'))
    row = best_move[1] + 1
    print(f"Best move is {col}{row} capturing {max_captures} stones")