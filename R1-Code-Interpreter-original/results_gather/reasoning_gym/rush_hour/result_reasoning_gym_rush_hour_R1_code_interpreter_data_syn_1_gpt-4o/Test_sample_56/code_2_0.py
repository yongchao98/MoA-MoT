# Initial board setup
board = [
    "..HBBB",
    "..HICC",
    "AAIJK",
    "GDDJK",
    "GEEJL",
    "FF...L"
]

# Function to move a car on the board
def move_car(board, car, direction, steps):
    new_board = [list(row) for row in board]
    if direction in ['right', 'left']:
        for i, row in enumerate(board):
            for j in range(len(row)):
                if row[j] == car:
                    if direction == 'right':
                        if j + steps < len(row) and all(new_board[i][k] == '.' for k in range(j + 1, j + steps + 1)):
                            new_board[i][j] = '.'
                            new_board[i][j + steps] = car
                    elif direction == 'left':
                        if j - steps >= 0 and all(new_board[i][k] == '.' for k in range(j - steps, j)):
                            new_board[i][j] = '.'
                            new_board[i][j - steps] = car
    elif direction in ['down', 'up']:
        for j in range(len(board[0])):
            for i in range(len(board)):
                if board[i][j] == car:
                    if direction == 'down':
                        if i + steps < len(board) and all(new_board[k][j] == '.' for k in range(i + 1, i + steps + 1)):
                            new_board[i][j] = '.'
                            new_board[i + steps][j] = car
                    elif direction == 'up':
                        if i - steps >= 0 and all(new_board[k][j] == '.' for k in range(i - steps, i)):
                            new_board[i][j] = '.'
                            new_board[i - steps][j] = car
    return [''.join(row) for row in new_board]

# Sequence of moves to solve the puzzle
moves = [
    ('J', 'down', 1),
    ('K', 'down', 1),
    ('L', 'down', 1),
    ('A', 'right', 2)
]

# Apply moves to the board
for car, direction, steps in moves:
    board = move_car(board, car, direction, steps)

# Print the final board state
for row in board:
    print(row)

# Output the moves in the required format
move_sequence = "J+1 K+1 L+1 A+2"
print(move_sequence)