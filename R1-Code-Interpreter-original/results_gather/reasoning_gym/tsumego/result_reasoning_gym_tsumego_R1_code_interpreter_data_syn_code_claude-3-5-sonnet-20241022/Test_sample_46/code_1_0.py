def print_board(board):
    print("   A B C D E F G H I J")
    for i in range(10):
        row = 10 - i
        print(f"{row:2d} {' '.join(board[i])}")

def is_valid(x, y):
    return 0 <= x < 10 and 0 <= y < 10

def get_neighbors(x, y):
    return [(x+dx, y+dy) for dx, dy in [(0,1), (0,-1), (1,0), (-1,0)] if is_valid(x+dx, y+dy)]

def count_captures(board, move_x, move_y):
    if board[move_y][move_x] != '.':
        return -1
    
    # Temporarily place the black stone
    board[move_y][move_x] = 'X'
    
    # Check for captures
    captured = set()
    checked = set()
    
    def check_group(x, y, stone):
        if (x, y) in checked:
            return False
        checked.add((x, y))
        
        if board[y][x] == '.':
            return False
        if board[y][x] != stone:
            return True
        
        group = {(x, y)}
        liberty = False
        stack = [(x, y)]
        
        while stack:
            cx, cy = stack.pop()
            for nx, ny in get_neighbors(cx, cy):
                if board[ny][nx] == '.':
                    liberty = True
                elif board[ny][nx] == stone and (nx, ny) not in group:
                    group.add((nx, ny))
                    stack.append((nx, ny))
        
        if not liberty:
            captured.update(group)
        return True

    # Check all adjacent positions
    for nx, ny in get_neighbors(move_x, move_y):
        if board[ny][nx] == 'O':
            check_group(nx, ny, 'O')
    
    # Remove the temporary stone
    board[move_y][move_x] = '.'
    
    return len(captured)

# Initialize the board
board = [['.'] * 10 for _ in range(10)]
# Set up the position
position = """
. . . . . . . . . .
. . . X O . . . . .
. . O . . . O . . .
. . X O . . . X . .
. . O . X . O . . .
. X O O O X O . . .
. . X O X . . . . .
. . . X . . . . . .
O . . . . . . . . .
. . . . . . . . . .
"""

# Parse the position
for i, row in enumerate(position.strip().split('\n')):
    for j, char in enumerate(row.split()):
        board[i][j] = char

# Find the best move
best_move = None
max_captures = 0

for y in range(10):
    for x in range(10):
        captures = count_captures(board, x, y)
        if captures > max_captures:
            max_captures = captures
            best_move = (x, y)

if best_move:
    col = chr(ord('A') + best_move[0])
    row = 10 - best_move[1]
    print(f"Best move is {col}{row} capturing {max_captures} stones")