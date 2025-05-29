def create_board():
    board = [['.'] * 13 for _ in range(13)]
    # Place White stones (O)
    white_positions = [
        (7, 12), (1, 10), (3, 9), (7, 6), (9, 6), (8, 4), (10, 5), (10, 2), (11, 2)
    ]
    # Place Black stones (X)
    black_positions = [
        (11, 9), (4, 8), (4, 7), (1, 5), (10, 6), (10, 4), (11, 6), (12, 5), (5, 1)
    ]
    
    for x, y in white_positions:
        board[y-1][x-1] = 'O'
    for x, y in black_positions:
        board[y-1][x-1] = 'X'
    
    return board

def count_liberties(board, x, y, visited=None):
    if visited is None:
        visited = set()
    
    if x < 0 or x >= 13 or y < 0 or y >= 13:
        return 0, set()
    
    if board[y][x] == '.':
        return 1, {(x, y)}
    
    if (x, y) in visited:
        return 0, set()
    
    if board[y][x] == '.':
        return 1, {(x, y)}
        
    visited.add((x, y))
    stone_type = board[y][x]
    liberties = set()
    
    for dx, dy in [(0, 1), (1, 0), (0, -1), (-1, 0)]:
        new_x, new_y = x + dx, y + dy
        if 0 <= new_x < 13 and 0 <= new_y < 13:
            if board[new_y][new_x] == '.':
                liberties.add((new_x, new_y))
            elif board[new_y][new_x] == stone_type:
                _, group_liberties = count_liberties(board, new_x, new_y, visited)
                liberties.update(group_liberties)
    
    return len(liberties), liberties

def find_best_capture():
    board = create_board()
    best_move = None
    max_captures = 0
    
    # Check each empty point
    for y in range(13):
        for x in range(13):
            if board[y][x] == '.':
                # Check if placing a stone here would capture any white groups
                captures = 0
                for dx, dy in [(0, 1), (1, 0), (0, -1), (-1, 0)]:
                    nx, ny = x + dx, y + dy
                    if 0 <= nx < 13 and 0 <= ny < 13 and board[ny][nx] == 'O':
                        liberties, _ = count_liberties(board, nx, ny)
                        if liberties == 1:  # This group would be captured
                            captures += 1
                if captures > max_captures:
                    max_captures = captures
                    # Convert to Go coordinates (A-M, 1-13)
                    best_move = (chr(x + 65), y + 1)

    return best_move

best_move = find_best_capture()
print(f"<<<{best_move[0]}{best_move[1]}>>>")