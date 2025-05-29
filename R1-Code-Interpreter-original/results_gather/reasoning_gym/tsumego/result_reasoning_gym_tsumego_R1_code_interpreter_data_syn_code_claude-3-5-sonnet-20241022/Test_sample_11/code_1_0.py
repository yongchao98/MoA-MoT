def print_board(board):
    print("   A B C D E F G H I")
    for i in range(9):
        row = 9 - i
        print(f" {row}", end=" ")
        for j in range(9):
            print(board[i][j], end=" ")
        print()

def get_liberties(board, x, y, color, visited=None):
    if visited is None:
        visited = set()
    
    if x < 0 or x >= 9 or y < 0 or y >= 9:
        return set()
    
    if board[x][y] == '.':
        return {(x, y)}
    
    if board[x][y] != color or (x, y) in visited:
        return set()
    
    visited.add((x, y))
    liberties = set()
    
    # Check all adjacent positions
    for dx, dy in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
        new_liberties = get_liberties(board, x + dx, y + dy, color, visited)
        liberties.update(new_liberties)
    
    return liberties

def find_groups(board, color):
    groups = []
    visited = set()
    
    for i in range(9):
        for j in range(9):
            if board[i][j] == color and (i, j) not in visited:
                group = set()
                stack = [(i, j)]
                while stack:
                    x, y = stack.pop()
                    if (x, y) not in visited and board[x][y] == color:
                        group.add((x, y))
                        visited.add((x, y))
                        for dx, dy in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
                            nx, ny = x + dx, y + dy
                            if 0 <= nx < 9 and 0 <= ny < 9:
                                stack.append((nx, ny))
                groups.append(group)
    return groups

# Initialize the board
board = [
    ['.', '.', 'O', '.', '.', '.', '.', '.', '.'],
    ['.', 'X', 'O', '.', '.', '.', '.', '.', '.'],
    ['X', '.', 'X', '.', '.', '.', '.', '.', '.'],
    ['O', 'O', 'O', 'X', '.', '.', '.', '.', '.'],
    ['X', 'O', 'O', '.', '.', '.', '.', '.', '.'],
    ['.', 'X', '.', '.', '.', '.', '.', '.', 'O'],
    ['.', 'X', '.', '.', '.', '.', 'X', '.', '.'],
    ['O', '.', 'O', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', 'O', '.', '.', '.', '.']
]

# Find all white groups and their liberties
white_groups = find_groups(board, 'O')
best_move = None
max_captures = 0

for group in white_groups:
    liberties = set()
    for x, y in group:
        liberties.update(get_liberties(board, x, y, 'O'))
    
    if len(liberties) == 1:
        liberty = liberties.pop()
        if len(group) > max_captures:
            max_captures = len(group)
            best_move = liberty

if best_move:
    # Convert to Go coordinates (A-I, 1-9)
    col = chr(best_move[1] + ord('A'))
    row = 9 - best_move[0]
    print(f"Best move: {col}{row} capturing {max_captures} stones")