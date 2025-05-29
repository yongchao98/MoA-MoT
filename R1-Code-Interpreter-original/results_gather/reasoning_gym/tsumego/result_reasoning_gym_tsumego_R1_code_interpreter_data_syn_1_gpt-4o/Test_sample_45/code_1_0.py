# Define the board
board = [
    ['.', '.', '.', '.', 'X', '.', '.', '.', '.', '.', '.'],
    ['.', 'O', '.', '.', 'X', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', 'X', '.'],
    ['.', '.', '.', 'X', '.', '.', '.', '.', 'O', 'X', 'X'],
    ['.', '.', '.', '.', '.', '.', '.', '.', 'X', 'O', 'O'],
    ['.', '.', '.', 'X', '.', '.', '.', '.', 'X', 'O', '.'],
    ['.', '.', 'X', '.', '.', '.', '.', '.', '.', 'X', '.'],
    ['.', '.', 'O', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', 'X', 'O'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.']
]

# Function to check if a move captures stones
def captures_stones(board, move, color):
    x, y = move
    opponent = 'O' if color == 'X' else 'X'
    directions = [(0, 1), (1, 0), (0, -1), (-1, 0)]
    captured = []

    for dx, dy in directions:
        nx, ny = x + dx, y + dy
        if 0 <= nx < len(board) and 0 <= ny < len(board[0]) and board[nx][ny] == opponent:
            # Check if this group is captured
            group = [(nx, ny)]
            visited = set(group)
            queue = list(group)
            is_captured = True

            while queue:
                cx, cy = queue.pop(0)
                for ddx, ddy in directions:
                    nnx, nny = cx + ddx, cy + ddy
                    if 0 <= nnx < len(board) and 0 <= nny < len(board[0]):
                        if board[nnx][nny] == '.':
                            is_captured = False
                        elif board[nnx][nny] == opponent and (nnx, nny) not in visited:
                            visited.add((nnx, nny))
                            queue.append((nnx, nny))
                            group.append((nnx, nny))

            if is_captured:
                captured.extend(group)

    return captured

# Check the move at H7
move = (4, 7)  # H7
captured_stones = captures_stones(board, move, 'X')
print(captured_stones)