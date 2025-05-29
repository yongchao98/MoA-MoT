def parse_board(board):
    white_groups = []
    visited = set()
    
    def dfs(x, y, group, liberties):
        if (x, y) in visited or board[y][x] != 'O':
            return
        visited.add((x, y))
        group.append((x, y))
        for dx, dy in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
            nx, ny = x + dx, y + dy
            if 0 <= nx < 10 and 0 <= ny < 10:
                if board[ny][nx] == '.':
                    liberties.add((nx, ny))
                elif board[ny][nx] == 'O':
                    dfs(nx, ny, group, liberties)
    
    for y in range(10):
        for x in range(10):
            if board[y][x] == 'O' and (x, y) not in visited:
                group = []
                liberties = set()
                dfs(x, y, group, liberties)
                white_groups.append((group, liberties))
    
    return white_groups

def simulate_move(board, move):
    x, y = move
    new_board = [list(row) for row in board]
    new_board[y][x] = 'X'
    return new_board

def count_captures(board, move):
    new_board = simulate_move(board, move)
    white_groups = parse_board(new_board)
    captures = 0
    for group, liberties in white_groups:
        if not liberties:
            captures += len(group)
    return captures

def find_best_move(board):
    best_move = None
    max_captures = 0
    for y in range(10):
        for x in range(10):
            if board[y][x] == '.':
                captures = count_captures(board, (x, y))
                if captures > max_captures:
                    max_captures = captures
                    best_move = (x, y)
    return best_move

board = [
    "..........",
    "..X.O.....",
    "..O...O...",
    "..XO...X..",
    "..O.X.O...",
    ".XOOOXO...",
    "..XOX.....",
    "...X......",
    "O.........",
    ".........."
]

best_move = find_best_move(board)
columns = "ABCDEFGHIJ"
if best_move:
    x, y = best_move
    print(f"<<<{columns[x]}{10-y}>>>")
else:
    print("No capturing move found.")