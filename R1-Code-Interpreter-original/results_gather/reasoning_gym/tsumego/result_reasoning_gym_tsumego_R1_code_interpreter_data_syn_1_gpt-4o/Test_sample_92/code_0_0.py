# Define the board
board = [
    ['.', '.', '.', '.', '.', '.', 'X', '.', 'X', '.', '.'],
    ['.', '.', '.', '.', '.', '.', 'X', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', 'X', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['X', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['X', '.', 'X', '.', 'O', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', 'X', '.', '.'],
    ['X', '.', '.', '.', '.', '.', 'X', 'O', 'O', '.', '.'],
    ['.', '.', '.', '.', '.', '.', 'X', 'O', 'O', '.', 'X'],
    ['.', '.', '.', '.', 'O', '.', 'X', 'O', 'X', '.', '.']
]

# Function to check if a move captures stones
def captures_stones(board, move, player):
    x, y = move
    opponent = 'O' if player == 'X' else 'X'
    captured = []

    # Check all four directions
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]
    for dx, dy in directions:
        nx, ny = x + dx, y + dy
        if 0 <= nx < len(board) and 0 <= ny < len(board[0]) and board[nx][ny] == opponent:
            # Check if the opponent's group is surrounded
            group = [(nx, ny)]
            queue = [(nx, ny)]
            surrounded = True
            while queue:
                cx, cy = queue.pop(0)
                for ddx, ddy in directions:
                    nnx, nny = cx + ddx, cy + ddy
                    if 0 <= nnx < len(board) and 0 <= nny < len(board[0]):
                        if board[nnx][nny] == '.':
                            surrounded = False
                        elif board[nnx][nny] == opponent and (nnx, nny) not in group:
                            group.append((nnx, nny))
                            queue.append((nnx, nny))
            if surrounded:
                captured.extend(group)

    return captured

# Check the move G2
move = (8, 6)  # G2
captured_stones = captures_stones(board, move, 'X')
print(len(captured_stones))