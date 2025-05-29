# Initial board setup
board = [
    ['.', '.', '.', 'F', '.', '.'],
    ['.', '.', 'x', 'F', 'G', '.'],
    ['A', 'A', 'E', 'F', 'G', '.'],
    ['.', '.', 'E', '.', 'H', '.'],
    ['.', 'D', 'B', 'B', 'H', '.'],
    ['.', 'D', 'C', 'C', 'x', '.']
]

# Function to print the board
def print_board(b):
    for row in b:
        print(''.join(row))
    print()

# Function to move a car
def move_car(b, car, direction, steps):
    # Find the car's position
    positions = [(i, j) for i in range(len(b)) for j in range(len(b[i])) if b[i][j] == car]
    if not positions:
        return b  # Car not found

    # Determine the orientation of the car
    if len(positions) < 2:
        return b  # Invalid car
    vertical = positions[0][0] != positions[1][0]

    # Move the car
    new_board = [row[:] for row in b]
    for i, j in positions:
        new_board[i][j] = '.'

    for i, j in positions:
        if vertical:
            new_i = i + steps if direction in ['down', 'up'] else i
            new_j = j
        else:
            new_i = i
            new_j = j + steps if direction in ['right', 'left'] else j
        new_board[new_i][new_j] = car

    return new_board

# Sequence of moves
moves = [
    ('G', 'down', 1),  # Move G down
    ('E', 'up', 1),    # Move E up
    ('A', 'right', 3)  # Move A to the exit
]

# Apply the moves
for car, direction, steps in moves:
    board = move_car(board, car, direction, steps)
    print_board(board)

# Final board state
print_board(board)