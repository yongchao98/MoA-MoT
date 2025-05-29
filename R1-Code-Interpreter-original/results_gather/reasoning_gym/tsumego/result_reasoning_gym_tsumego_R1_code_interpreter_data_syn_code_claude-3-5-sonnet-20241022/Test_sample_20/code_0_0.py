def print_board(board):
    print("   A B C D E F G H I J K")
    for i in range(11, 0, -1):
        row = f"{i:2d} "
        for j in range(11):
            row += board[i-1][j] + " "
        print(row)

def is_valid_move(x, y):
    return 0 <= x < 11 and 0 <= y < 11

def get_liberties(board, x, y, color, visited):
    if not is_valid_move(x, y):
        return 0
    if (x, y) in visited:
        return 0
    if board[y][x] == '.':
        return 1
    if board[y][x] != color:
        return 0
    
    visited.add((x, y))
    liberties = 0
    for dx, dy in [(0,1), (0,-1), (1,0), (-1,0)]:
        liberties += get_liberties(board, x+dx, y+dy, color, visited)
    return liberties

def analyze_position():
    board = [
        ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
        ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
        ['.', '.', '.', '.', '.', '.', 'X', '.', '.', 'X', '.'],
        ['.', '.', '.', 'O', '.', '.', '.', 'X', 'X', 'O', 'X'],
        ['.', '.', 'O', '.', '.', '.', '.', 'X', 'O', 'O', '.'],
        ['.', '.', '.', '.', '.', '.', 'X', '.', 'X', 'O', 'X'],
        ['.', '.', '.', '.', '.', '.', '.', '.', '.', 'X', '.'],
        ['.', 'X', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
        ['.', '.', 'O', '.', '.', '.', '.', '.', '.', '.', '.'],
        ['.', '.', '.', '.', '.', '.', '.', '.', '.', 'X', '.'],
        ['.', '.', '.', '.', '.', '.', '.', '.', 'X', 'O', '.']
    ]

    # Check key positions for potential moves
    critical_moves = []
    for y in range(11):
        for x in range(11):
            if board[y][x] == '.':
                # Try move
                board[y][x] = 'X'
                # Check if any adjacent white groups are in atari
                for dx, dy in [(0,1), (0,-1), (1,0), (-1,0)]:
                    nx, ny = x+dx, y+dy
                    if is_valid_move(nx, ny) and board[ny][nx] == 'O':
                        liberties = get_liberties(board, nx, ny, 'O', set())
                        if liberties <= 2:  # Looking for moves that put white in atari or capture
                            col = chr(ord('A') + x)
                            row = 11 - y
                            critical_moves.append(f"{col}{row}")
                board[y][x] = '.'

    print("Critical moves to consider:", critical_moves)

analyze_position()