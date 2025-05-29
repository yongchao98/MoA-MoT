# Initial board setup
board = [
    ['x', '.', 'G', 'B', 'B', 'J'],
    ['C', 'C', 'G', 'H', 'I', 'J'],
    ['F', 'A', 'A', 'H', 'I', 'K'],
    ['F', 'D', 'D', '.', 'I', 'K'],
    ['E', 'E', '.', 'x', '.', '.'],
    ['.', '.', '.', '.', '.', '.']
]

# Function to move a car
def move_car(board, car, direction, steps):
    if direction == 'right':
        for _ in range(steps):
            for row in board:
                for i in range(len(row) - 1, -1, -1):
                    if row[i] == car and i + 1 < len(row) and row[i + 1] == '.':
                        row[i], row[i + 1] = row[i + 1], row[i]
    elif direction == 'down':
        for _ in range(steps):
            for i in range(len(board) - 1, -1, -1):
                for j in range(len(board[i])):
                    if board[i][j] == car and i + 1 < len(board) and board[i + 1][j] == '.':
                        board[i][j], board[i + 1][j] = board[i + 1][j], board[i]

# Move sequence
move_car(board, 'H', 'down', 1)
move_car(board, 'I', 'down', 1)
move_car(board, 'J', 'down', 1)
move_car(board, 'A', 'right', 3)

# Print the final board state
for row in board:
    print(''.join(str(cell) for cell in row))