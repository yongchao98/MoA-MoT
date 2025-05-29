# Define the board
board = [
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['O', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['O', '.', '.', '.', '.', '.', 'X', 'X', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', 'X', '.', '.'],
    ['X', '.', '.', '.', '.', '.', '.', '.', 'X', 'O', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', 'X', 'X', 'O', 'O', 'X'],
    ['.', '.', 'O', '.', '.', 'O', '.', 'O', '.', 'X', 'X', '.'],
    ['.', 'X', '.', '.', '.', '.', '.', 'X', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', 'X']
]

# Function to check if a group of stones is captured
def is_captured(board, x, y, color):
    visited = set()
    stack = [(x, y)]
    liberties = 0

    while stack:
        cx, cy = stack.pop()
        if (cx, cy) in visited:
            continue
        visited.add((cx, cy))

        # Check all adjacent positions
        for dx, dy in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
            nx, ny = cx + dx, cy + dy
            if 0 <= nx < len(board) and 0 <= ny < len(board[0]):
                if board[nx][ny] == '.':
                    liberties += 1
                elif board[nx][ny] == color and (nx, ny) not in visited:
                    stack.append((nx, ny))

    return liberties == 0

# Check all white stones to see if they can be captured
for i in range(len(board)):
    for j in range(len(board[0])):
        if board[i][j] == 'O':
            # Check if placing a black stone around this group captures it
            for dx, dy in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
                ni, nj = i + dx, j + dy
                if 0 <= ni < len(board) and 0 <= nj < len(board[0]) and board[ni][nj] == '.':
                    # Temporarily place a black stone
                    board[ni][nj] = 'X'
                    if is_captured(board, i, j, 'O'):
                        print(f"Placing at {ni+1}, {chr(nj+65)} captures a group.")
                    # Remove the temporary black stone
                    board[ni][nj] = '.'