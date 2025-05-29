# Define the board
board = [
    ['.', '.', '.', '.', 'O', '.', '.', 'O', '.'],
    ['X', '.', '.', '.', 'X', 'O', 'X', '.', 'X'],
    ['.', '.', 'X', 'X', 'O', 'O', 'X', '.', '.'],
    ['.', '.', '.', 'X', 'O', '.', 'X', '.', '.'],
    ['.', '.', '.', '.', 'X', '.', 'O', '.', '.'],
    ['.', '.', '.', '.', 'O', '.', 'O', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.']
]

# Function to find a group of connected stones and their liberties
def find_group(board, start, color):
    stack = [start]
    visited = set()
    group = []
    liberties = 0

    while stack:
        x, y = stack.pop()
        if (x, y) in visited:
            continue
        visited.add((x, y))
        group.append((x, y))

        for dx, dy in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
            nx, ny = x + dx, y + dy
            if 0 <= nx < 9 and 0 <= ny < 9:
                if board[nx][ny] == '.':
                    liberties += 1
                elif board[nx][ny] == color and (nx, ny) not in visited:
                    stack.append((nx, ny))

    return group, liberties

# Function to check if a move captures any stones
def check_capture(board, move, color):
    x, y = move
    opponent = 'O' if color == 'X' else 'X'
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]
    captured = []

    for dx, dy in directions:
        nx, ny = x + dx, y + dy
        if 0 <= nx < 9 and 0 <= ny < 9 and board[nx][ny] == opponent:
            # Check if the opponent's group is captured
            group, liberties = find_group(board, (nx, ny), opponent)
            if liberties == 0:
                captured.extend(group)

    return captured

# Check the move G8
move = (1, 6)  # G8 in 0-indexed coordinates
captured_stones = check_capture(board, move, 'X')
print(captured_stones)