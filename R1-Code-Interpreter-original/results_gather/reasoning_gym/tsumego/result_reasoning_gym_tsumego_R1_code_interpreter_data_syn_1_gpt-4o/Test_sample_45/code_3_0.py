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

    def has_liberty(x, y, visited):
        if (x, y) in visited:
            return False
        visited.add((x, y))
        for dx, dy in directions:
            nx, ny = x + dx, y + dy
            if 0 <= nx < len(board) and 0 <= ny < len(board[0]):
                if board[nx][ny] == '.':
                    return True
                if board[nx][ny] == opponent and has_liberty(nx, ny, visited):
                    return True
        return False

    for dx, dy in directions:
        nx, ny = x + dx, y + dy
        if 0 <= nx < len(board) and 0 <= ny < len(board[0]) and board[nx][ny] == opponent:
            visited = set()
            if not has_liberty(nx, ny, visited):
                captured.extend(visited)

    return captured

# Check the move at H7
move = (6, 7)  # H7
captured_stones = captures_stones(board, move, 'X')
print(captured_stones)