def is_valid(x, y):
    return 0 <= x < 9 and 0 <= y < 9

def get_liberties(board, x, y, color, visited):
    if not is_valid(x, y) or (x, y) in visited:
        return set()
    if board[y][x] == '.':
        return {(x, y)}
    if board[y][x] != color:
        return set()
    
    visited.add((x, y))
    liberties = set()
    for dx, dy in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
        new_x, new_y = x + dx, y + dy
        liberties.update(get_liberties(board, new_x, new_y, color, visited))
    return liberties

def analyze_position():
    board = [
        ['.', '.', '.', '.', '.', '.', '.', '.', '.'],
        ['.', '.', '.', '.', '.', '.', '.', 'O', '.'],
        ['.', '.', 'X', '.', '.', '.', 'O', '.', '.'],
        ['.', '.', '.', '.', '.', '.', 'O', '.', '.'],
        ['.', '.', '.', 'X', '.', '.', 'O', '.', '.'],
        ['.', 'O', '.', '.', '.', '.', 'O', 'X', '.'],
        ['X', '.', 'X', '.', '.', '.', '.', '.', '.'],
        ['O', 'O', 'O', 'X', '.', '.', '.', 'X', '.'],
        ['X', 'O', 'X', 'O', '.', '.', 'X', '.', '.']
    ]
    
    # Check each empty point for potential capturing moves
    for y in range(9):
        for x in range(9):
            if board[y][x] == '.':
                # Try placing a black stone here
                board[y][x] = 'X'
                # Check surrounding white groups
                for dx, dy in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
                    nx, ny = x + dx, y + dy
                    if is_valid(nx, ny) and board[ny][nx] == 'O':
                        liberties = get_liberties(board, nx, ny, 'O', set())
                        if not liberties:
                            print(f"Capturing move found at {chr(65+x)}{9-y}")
                board[y][x] = '.'

analyze_position()