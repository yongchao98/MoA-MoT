def is_valid(x, y):
    return 0 <= x < 11 and 0 <= y < 11

def count_liberties(board, x, y, visited):
    if not is_valid(x, y) or (x, y) in visited:
        return 0
    if board[y][x] == '.':
        return 1
    if board[y][x] != 'O':
        return 0
    
    visited.add((x, y))
    liberties = 0
    for dx, dy in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
        liberties += count_liberties(board, x + dx, y + dy, visited)
    return liberties

def analyze_position():
    board = [
        ['.', '.', '.', '.', '.', '.', '.', '.', '.', 'O', '.'],
        ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', 'O'],
        ['.', '.', 'X', '.', '.', 'X', 'X', '.', '.', '.', '.'],
        ['.', '.', '.', '.', 'X', 'O', 'X', 'X', '.', '.', '.'],
        ['.', '.', 'X', 'X', 'O', 'O', '.', '.', '.', 'X', '.'],
        ['O', '.', '.', '.', 'X', 'O', 'X', '.', '.', '.', '.'],
        ['.', '.', '.', '.', '.', 'X', '.', '.', '.', '.', '.'],
        ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
        ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
        ['.', '.', '.', '.', '.', '.', '.', '.', 'O', '.', '.'],
        ['X', '.', '.', '.', '.', '.', '.', '.', '.', '.', 'X']
    ]
    
    best_move = None
    max_captures = 0
    
    # Check each empty position
    for y in range(11):
        for x in range(11):
            if board[y][x] == '.':
                # Try placing a black stone here
                board[y][x] = 'X'
                
                # Check surrounding white groups
                total_captures = 0
                for dx, dy in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
                    nx, ny = x + dx, y + dy
                    if is_valid(nx, ny) and board[ny][nx] == 'O':
                        liberties = count_liberties(board, nx, ny, set())
                        if liberties == 0:
                            # Count stones in this captured group
                            visited = set()
                            def count_group(x, y):
                                if not is_valid(x, y) or (x, y) in visited or board[y][x] != 'O':
                                    return 0
                                visited.add((x, y))
                                count = 1
                                for dx, dy in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
                                    count += count_group(x + dx, y + dy)
                                return count
                            total_captures += count_group(nx, ny)
                
                if total_captures > max_captures:
                    max_captures = total_captures
                    best_move = (x, y)
                
                # Reset the position
                board[y][x] = '.'
    
    # Convert to board coordinates
    if best_move:
        col = chr(ord('A') + best_move[0])
        row = 11 - best_move[1]
        print(f"{col}{row} (captures {max_captures} stones)")

analyze_position()