# Define the board
board = [
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', 'X', '.', 'O', '.', '.', '.', '.'],
    ['X', 'O', 'O', '.', '.', '.', '.', '.', '.', '.'],
    ['.', 'X', 'X', '.', '.', '.', '.', '.', '.', '.'],
    ['O', 'O', 'O', 'X', '.', '.', 'O', '.', '.', '.'],
    ['X', 'O', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', 'X', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', 'X', '.', '.', '.', '.', '.', '.', 'X', '.'],
    ['.', '.', 'X', '.', '.', '.', '.', 'X', '.', '.'],
    ['.', 'X', '.', '.', '.', '.', '.', '.', '.', '.']
]

# Function to check if a position is within the board
def is_within_board(x, y):
    return 0 <= x < 10 and 0 <= y < 10

# Function to get the liberties of a group
def get_liberties(x, y, color):
    visited = set()
    stack = [(x, y)]
    liberties = set()
    while stack:
        cx, cy = stack.pop()
        if (cx, cy) in visited:
            continue
        visited.add((cx, cy))
        for dx, dy in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
            nx, ny = cx + dx, cy + dy
            if is_within_board(nx, ny):
                if board[ny][nx] == '.':
                    liberties.add((nx, ny))
                elif board[ny][nx] == color and (nx, ny) not in visited:
                    stack.append((nx, ny))
    return liberties

# Function to find the best move for Black
def find_best_move():
    # Focus on the group at A6, B6, C6
    target_group = [(0, 4), (1, 4), (2, 4)]
    liberties = set()
    for x, y in target_group:
        liberties.update(get_liberties(x, y, 'O'))
    
    # Check if placing a black stone at any liberty captures the group
    for lx, ly in liberties:
        captures = sum(1 for dx, dy in [(-1, 0), (1, 0), (0, -1), (0, 1)]
                       if is_within_board(lx + dx, ly + dy) and board[ly + dy][lx + dx] == 'O')
        if captures == len(target_group):
            return (lx, ly)
    return None

# Find and print the best move
best_move = find_best_move()
if best_move:
    column = chr(best_move[0] + ord('A'))
    row = 10 - best_move[1]
    print(f"<<<{column}{row}>>>")
else:
    print("No capturing move found")