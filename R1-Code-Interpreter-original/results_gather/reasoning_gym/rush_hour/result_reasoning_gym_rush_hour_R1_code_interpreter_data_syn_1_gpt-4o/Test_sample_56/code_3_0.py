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
    car_positions = [(i, j) for i, row in enumerate(board) for j, cell in enumerate(row) if cell == car]
    
    if not car_positions:
        return board  # Car not found, return the board as is
    
    if direction in ['right', 'left']:
        row = car_positions[0][0]
        if direction == 'right':
            max_j = max(j for _, j in car_positions)
            if max_j + steps < len(board[0]) and all(new_board[row][j] == '.' for j in range(max_j + 1, max_j + steps + 1)):
                for _, j in car_positions:
                    new_board[row][j] = '.'
                for j in range(max_j + 1, max_j + steps + 1):
                    new_board[row][j] = car
        elif direction == 'left':
            min_j = min(j for _, j in car_positions)
            if min_j - steps >= 0 and all(new_board[row][j] == '.' for j in range(min_j - steps, min_j)):
                for _, j in car_positions:
                    new_board[row][j] = '.'
                for j in range(min_j - steps, min_j):
                    new_board[row][j] = car
    elif direction in ['down', 'up']:
        col = car_positions[0][1]
        if direction == 'down':
            max_i = max(i for i, _ in car_positions)
            if max_i + steps < len(board) and all(new_board[i][col] == '.' for i in range(max_i + 1, max_i + steps + 1)):
                for i, _ in car_positions:
                    new_board[i][col] = '.'
                for i in range(max_i + 1, max_i + steps + 1):
                    new_board[i][col] = car
        elif direction == 'up':
            min_i = min(i for i, _ in car_positions)
            if min_i - steps >= 0 and all(new_board[i][col] == '.' for i in range(min_i - steps, min_i)):
                for i, _ in car_positions:
                    new_board[i][col] = '.'
                for i in range(min_i - steps, min_i):
                    new_board[i][col] = car
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