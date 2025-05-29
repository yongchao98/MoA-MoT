def simulate_move(board, move, player):
    # Simulate the move on the board
    x, y = move
    board[y][x] = player
    return board

def count_liberties(board, x, y):
    # Count liberties for a group of stones
    visited = set()
    stack = [(x, y)]
    liberties = 0
    while stack:
        cx, cy = stack.pop()
        if (cx, cy) in visited:
            continue
        visited.add((cx, cy))
        for nx, ny in [(cx-1, cy), (cx+1, cy), (cx, cy-1), (cx, cy+1)]:
            if 0 <= nx < len(board[0]) and 0 <= ny < len(board):
                if board[ny][nx] == '.':
                    liberties += 1
                elif board[ny][nx] == board[cy][cx]:
                    stack.append((nx, ny))
    return liberties

def find_best_move(board):
    # Check potential moves and their impact
    moves = [(3, 8), (4, 8)]  # D8, E8
    best_move = None
    max_captured = 0
    for move in moves:
        new_board = [row[:] for row in board]
        simulate_move(new_board, move, 'X')
        captured = 0
        for x, y in [(2, 8), (3, 7), (3, 9)]:
            if new_board[y][x] == 'O' and count_liberties(new_board, x, y) == 0:
                captured += 1
        if captured > max_captured:
            max_captured = captured
            best_move = move
    return best_move

# Initial board setup
board = [
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', 'O', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', 'X', '.', '.', 'X', '.', '.', '.', '.', '.', '.'],
    ['.', 'X', 'O', 'X', '.', '.', '.', '.', '.', 'O', '.', '.'],
    ['X', 'O', 'O', '.', '.', '.', '.', 'O', '.', '.', '.', '.'],
    ['.', 'X', 'O', 'X', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', 'X', '.', '.', '.', '.', '.', '.', '.', 'X', 'X'],
    ['.', '.', '.', '.', 'O', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['X', '.', '.', '.', '.', '.', '.', '.', '.', '.', 'X', '.'],
    ['.', '.', '.', 'X', 'O', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.']
]

best_move = find_best_move(board)
print(best_move)