def create_board():
    board = [['.'] * 10 for _ in range(10)]
    # Place black stones (X)
    black_stones = [(9,2), (8,0), (7,1), (7,3), (7,9), (6,2), (6,9), (5,9), (4,9), (1,1), (1,8)]
    # Place white stones (O)
    white_stones = [(8,1), (8,2), (7,2), (5,5), (5,8), (3,3), (2,8)]
    
    for x, y in black_stones:
        board[x][y] = 'X'
    for x, y in white_stones:
        board[x][y] = 'O'
    return board

def check_capture_sequence(board, x, y):
    # Check if playing at (x,y) leads to immediate or subsequent captures
    stones_captured = 0
    groups_to_check = []
    
    # Check adjacent positions for white stones
    for dx, dy in [(0,1), (0,-1), (1,0), (-1,0)]:
        nx, ny = x + dx, y + dy
        if 0 <= nx < 10 and 0 <= ny < 10 and board[nx][ny] == 'O':
            liberties = 0
            group_size = 0
            visited = set()
            stack = [(nx, ny)]
            
            while stack:
                cx, cy = stack.pop()
                if (cx, cy) in visited:
                    continue
                visited.add((cx, cy))
                group_size += 1
                
                for ddx, ddy in [(0,1), (0,-1), (1,0), (-1,0)]:
                    nnx, nny = cx + ddx, cy + ddy
                    if 0 <= nnx < 10 and 0 <= nny < 10:
                        if board[nnx][nny] == '.':
                            if (nnx, nny) != (x, y):  # Don't count the move we're checking
                                liberties += 1
                        elif board[nnx][nny] == 'O' and (nnx, nny) not in visited:
                            stack.append((nnx, nny))
            
            if liberties == 0:  # This group would be captured
                stones_captured += group_size

    print(f"Move at {chr(y + ord('A'))}{10-x} captures {stones_captured} stones")
    return stones_captured

def find_best_capture():
    board = create_board()
    best_move = None
    max_captures = 0
    
    # Try all empty positions
    for i in range(10):
        for j in range(10):
            if board[i][j] == '.':
                captures = check_capture_sequence(board, i, j)
                if captures > max_captures:
                    max_captures = captures
                    best_move = (i, j)
    
    if best_move:
        print(f"\nBest move: {chr(best_move[1] + ord('A'))}{10-best_move[0]} capturing {max_captures} stones")

find_best_capture()