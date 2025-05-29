def create_board():
    board = [['.'] * 12 for _ in range(12)]
    # Place black stones (X)
    black_stones = [(11,0), (10,4), (10,6), (10,7), (9,5), (9,8), (9,11),
                    (8,1), (8,5), (7,1), (7,2), (7,6), (5,5), (1,10)]
    # Place white stones (O)
    white_stones = [(12,8), (9,4), (9,6), (9,7), (8,6), (5,0), (1,7)]
    
    for x, y in black_stones:
        board[12-x][y] = 'X'
    for x, y in white_stones:
        board[12-x][y] = 'O'
    return board

def is_valid(x, y):
    return 0 <= x < 12 and 0 <= y < 12

def count_captures(board, x, y):
    if board[x][y] != '.':
        return 0
    
    directions = [(0,1), (1,0), (0,-1), (-1,0)]
    captures = 0
    
    # Try placing a black stone and check for captures
    board[x][y] = 'X'
    
    def is_captured(x, y, visited):
        if not is_valid(x, y) or board[x][y] == 'X':
            return False
        if board[x][y] == '.':
            return False
        if (x,y) in visited:
            return True
            
        visited.add((x,y))
        for dx, dy in directions:
            nx, ny = x + dx, y + dy
            if not is_captured(nx, ny, visited):
                return False
        return True
    
    # Check each direction for potential captures
    for dx, dy in directions:
        nx, ny = x + dx, y + dy
        if is_valid(nx, ny) and board[nx][ny] == 'O':
            if is_captured(nx, ny, set()):
                # Count connected white stones
                def count_group(x, y, counted):
                    if not is_valid(x, y) or board[x][y] != 'O' or (x,y) in counted:
                        return 0
                    counted.add((x,y))
                    count = 1
                    for dx, dy in directions:
                        count += count_group(x+dx, y+dy, counted)
                    return count
                captures += count_group(nx, ny, set())
    
    board[x][y] = '.'
    return captures

def find_best_move():
    board = create_board()
    best_move = None
    max_captures = 0
    
    for i in range(12):
        for j in range(12):
            if board[i][j] == '.':
                captures = count_captures(board, i, j)
                if captures > max_captures:
                    max_captures = captures
                    best_move = (i, j)
    
    if best_move:
        # Convert to Go coordinates (A-L, 1-12)
        col = chr(best_move[1] + ord('A'))
        row = 12 - best_move[0]
        return f"{col}{row}", max_captures
    return None, 0

move, captures = find_best_move()
print(f"Best move: {move}")
print(f"Captures: {captures}")