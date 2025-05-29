def create_board():
    board = [['.'] * 11 for _ in range(11)]
    # Place black stones (X)
    black_positions = [(10,'A'), (9,'B'), (8,'A'), (8,'C'), (7,'D'), 
                      (6,'A'), (6,'C'), (6,'J'), (5,'B'), (5,'K'),
                      (4,'D'), (4,'I'), (2,'F')]
    # Place white stones (O)
    white_positions = [(10,'G'), (8,'B'), (8,'I'), (7,'A'), (7,'B'),
                      (6,'B'), (6,'F'), (5,'G'), (3,'J')]
    
    # Convert to 0-based indexing and place stones
    for row, col in black_positions:
        board[11-row][ord(col)-ord('A')] = 'X'
    for row, col in white_positions:
        board[11-row][ord(col)-ord('A')] = 'O'
    return board

def get_neighbors(x, y, size=11):
    return [(x+dx, y+dy) for dx, dy in [(0,1), (0,-1), (1,0), (-1,0)]
            if 0 <= x+dx < size and 0 <= y+dy < size]

def find_group(board, x, y, stone_type, visited=None):
    if visited is None:
        visited = set()
    
    if (x, y) in visited or board[x][y] != stone_type:
        return set()
    
    group = {(x, y)}
    visited.add((x, y))
    
    for nx, ny in get_neighbors(x, y):
        if board[nx][ny] == stone_type:
            group.update(find_group(board, nx, ny, stone_type, visited))
    
    return group

def count_liberties(board, group):
    liberties = set()
    for x, y in group:
        for nx, ny in get_neighbors(x, y):
            if board[nx][ny] == '.':
                liberties.add((nx, ny))
    return liberties

def evaluate_move(board, x, y):
    if board[x][y] != '.':
        return 0
    
    # Try the move
    board[x][y] = 'X'
    captured = 0
    
    # Check all neighboring white groups
    checked = set()
    for nx, ny in get_neighbors(x, y):
        if board[nx][ny] == 'O' and (nx, ny) not in checked:
            group = find_group(board, nx, ny, 'O')
            checked.update(group)
            if not count_liberties(board, group):
                captured += len(group)
    
    # Restore board
    board[x][y] = '.'
    return captured

def find_best_move(board):
    best_score = 0
    best_move = None
    
    for i in range(11):
        for j in range(11):
            if board[i][j] == '.':
                score = evaluate_move(board, i, j)
                if score > best_score:
                    best_score = score
                    best_move = (i, j)
    
    if best_move:
        return chr(best_move[1] + ord('A')) + str(11 - best_move[0])
    return None

board = create_board()
best_move = find_best_move(board)
print(best_move)