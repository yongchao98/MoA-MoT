# Initial board setup
board = [
    ['.', '.', 'D', '.', '.', '.'],
    ['.', '.', 'D', 'B', 'B', 'B'],
    ['.', 'A', 'A', 'E', '.', 'F'],
    ['.', '.', '.', 'E', '.', 'F'],
    ['.', 'C', 'C', 'x', '.', 'F'],
    ['.', '.', '.', '.', '.', '.']
]

# Function to move a car
def move_car(board, car, direction, steps):
    for _ in range(steps):
        if direction == 'right':
            for i in range(len(board)):
                for j in range(len(board[i]) - 1, -1, -1):
                    if board[i][j] == car and j + 1 < len(board[i]) and board[i][j + 1] == '.':
                        board[i][j], board[i][j + 1] = board[i][j + 1], board[i][j]
        elif direction == 'left':
            for i in range(len(board)):
                for j in range(len(board[i])):
                    if board[i][j] == car and j - 1 >= 0 and board[i][j - 1] == '.':
                        board[i][j], board[i][j - 1] = board[i][j - 1], board[i][j]
        elif direction == 'down':
            for i in range(len(board) - 1, -1, -1):
                for j in range(len(board[i])):
                    if board[i][j] == car and i + 1 < len(board) and board[i + 1][j] == '.':
                        board[i][j], board[i + 1][j] = board[i + 1][j], board[i][j]
        elif direction == 'up':
            for i in range(len(board)):
                for j in range(len(board[i])):
                    if board[i][j] == car and i - 1 >= 0 and board[i - 1][j] == '.':
                        board[i][j], board[i - 1][j] = board[i - 1][j], board[i][j]

# Sequence of moves
move_car(board, 'F', 'down', 1)
move_car(board, 'E', 'right', 1)
move_car(board, 'A', 'right', 3)

# Print the final board state
for row in board:
    print(''.join(row))