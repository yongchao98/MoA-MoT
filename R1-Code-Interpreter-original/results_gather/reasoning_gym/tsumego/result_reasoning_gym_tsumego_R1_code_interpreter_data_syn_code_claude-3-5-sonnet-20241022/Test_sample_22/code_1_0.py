def create_board():
    board = [['.'] * 10 for _ in range(10)]
    # Place stones according to the given position
    black = [(7,1), (8,2), (7,3), (7,2), (5,0), (5,2), (4,0), (4,1), (3,6), (3,8)]
    white = [(6,0), (6,1), (6,2), (5,1), (4,4), (2,5), (2,7), (1,3), (1,4), (0,9)]
    
    for x, y in black:
        board[x][y] = 'X'
    for x, y in white:
        board[x][y] = 'O'
    return board

def is_valid_pos(x, y):
    return 0 <= x < 10 and 0 <= y < 10

def get_group(board, x, y, color, visited=None):
    if visited is None:
        visited = set()
    if not is_valid_pos(x, y) or board[x][y] != color or (x,y) in visited:
        return visited
    
    visited.add((x,y))
    for dx, dy in [(0,1), (0,-1), (1,0), (-1,0)]:
        new_x, new_y = x + dx, y + dy
        get_group(board, new_x, new_y, color, visited)
    return visited

def count_liberties_group(board, group):
    liberties = set()
    for x, y in group:
        for dx, dy in [(0,1), (0,-1), (1,0), (-1,0)]:
            new_x, new_y = x + dx, y + dy
            if is_valid_pos(new_x, new_y) and board[new_x][new_y] == '.':
                liberties.add((new_x, new_y))
    return len(liberties)

def evaluate_move(board, x, y):
    if board[x][y] != '.':
        return 0
    
    # Try the move
    board[x][y] = 'X'
    captured = 0
    
    # Check all adjacent positions for white groups
    for dx, dy in [(0,1), (0,-1), (1,0), (-1,0)]:
        new_x, new_y = x + dx, y + dy
        if is_valid_pos(new_x, new_y) and board[new_x][new_y] == 'O':
            group = get_group(board, new_x, new_y, 'O')
            if count_liberties_group(board, group) == 0:
                captured += len(group)
    
    # Undo the move
    board[x][y] = '.'
    return captured

board = create_board()
best_move = None
max_captured = 0

for i in range(10):
    for j in range(10):
        if board[i][j] == '.':
            captured = evaluate_move(board, i, j)
            if captured > max_captured:
                max_captured = captured
                best_move = (i, j)

if best_move:
    print(f"Best move: ({chr(65+best_move[1])}{10-best_move[0]}) capturing {max_captured} stones")