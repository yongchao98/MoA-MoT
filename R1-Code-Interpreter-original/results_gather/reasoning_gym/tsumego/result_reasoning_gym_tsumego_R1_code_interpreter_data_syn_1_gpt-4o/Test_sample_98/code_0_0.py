def count_captured_stones(board, move):
    # Simulate the move and count captured stones
    def is_captured(x, y, color):
        # Check if a group of stones is captured
        visited = set()
        stack = [(x, y)]
        liberties = 0
        while stack:
            cx, cy = stack.pop()
            if (cx, cy) in visited:
                continue
            visited.add((cx, cy))
            for nx, ny in [(cx-1, cy), (cx+1, cy), (cx, cy-1), (cx, cy+1)]:
                if 0 <= nx < len(board) and 0 <= ny < len(board[0]):
                    if board[nx][ny] == '.':
                        liberties += 1
                    elif board[nx][ny] == color:
                        stack.append((nx, ny))
        return liberties == 0

    x, y = move
    board[x][y] = 'X'  # Place the black stone
    captured = 0
    for nx, ny in [(x-1, y), (x+1, y), (x, y-1), (x, y+1)]:
        if 0 <= nx < len(board) and 0 <= ny < len(board[0]) and board[nx][ny] == 'O':
            if is_captured(nx, ny, 'O'):
                captured += 1
    return captured

# Board representation
board = [
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['X', '.', 'O', '.', '.', '.', '.', 'O', '.', '.', 'O', '.', '.'],
    ['.', 'X', '.', 'X', 'O', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['X', 'O', 'O', 'O', 'X', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', 'X', 'O', 'X', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', 'X', '.', '.', 'X', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', 'O', '.', '.', 'O', 'X', 'X'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', 'X', '.', '.', '.', '.', '.', '.', '.']
]

# Evaluate potential moves
move_f7 = (7, 5)  # F7
move_b8 = (6, 1)  # B8

captured_f7 = count_captured_stones(board, move_f7)
captured_b8 = count_captured_stones(board, move_b8)

print("Captured stones at F7:", captured_f7)
print("Captured stones at B8:", captured_b8)