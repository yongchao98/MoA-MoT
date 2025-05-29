def create_board():
    board = [['.'] * 13 for _ in range(13)]
    # Place White stones (O)
    white_positions = [(0,0), (3,1), (2,1), (3,2), (4,2), (5,2), (3,3), (4,3), 
                      (3,4), (11,5), (9,13)]
    # Place Black stones (X)
    black_positions = [(2,2), (2,3), (5,12), (2,11), (7,11), (9,8), (9,7), (4,1)]
    
    for x, y in white_positions:
        board[y][x] = 'O'
    for x, y in black_positions:
        board[y][x] = 'X'
    return board

def is_valid(x, y):
    return 0 <= x < 13 and 0 <= y < 13

def count_captures(board, x, y):
    if board[y][x] != '.':
        return 0
    
    directions = [(0,1), (1,0), (0,-1), (-1,0)]
    captured = set()
    
    for dx, dy in directions:
        stones = []
        nx, ny = x + dx, y + dy
        
        while is_valid(nx, ny) and board[ny][nx] == 'O':
            stones.append((nx, ny))
            nx += dx
            ny += dy
            
        if stones and is_valid(nx, ny) and board[ny][nx] == 'X':
            captured.update(stones)
    
    return len(captured)

def find_best_move():
    board = create_board()
    best_move = None
    max_captures = 0
    
    for y in range(13):
        for x in range(13):
            if board[y][x] == '.':
                captures = count_captures(board, x, y)
                if captures > max_captures:
                    max_captures = captures
                    best_move = (x, y)
    
    if best_move:
        col = chr(ord('A') + best_move[0])
        row = 13 - best_move[1]
        print(f"Best move: {col}{row} with {max_captures} captures")

find_best_move()