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
    for i, row in enumerate(board):
        for j, cell in enumerate(row):
            if cell == car:
                new_board[i][j] = '.'
                if direction == 'right':
                    new_board[i][j + steps] = car
                elif direction == 'down':
                    new_board[i + steps][j] = car
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