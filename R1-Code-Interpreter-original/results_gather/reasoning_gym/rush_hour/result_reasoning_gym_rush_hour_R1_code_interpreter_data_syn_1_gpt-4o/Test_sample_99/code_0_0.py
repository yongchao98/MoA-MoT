# Initial board setup
board = [
    ['B', 'B', 'C', 'C', '.', 'x'],
    ['.', '.', '.', 'K', 'D', 'D'],
    ['x', 'A', 'A', 'K', '.', 'L'],
    ['E', 'E', 'J', 'F', 'F', 'L'],
    ['.', '.', 'J', '.', 'G', 'G'],
    ['.', 'H', 'H', 'H', 'I', 'I']
]

# Function to print the board
def print_board(b):
    for row in b:
        print(''.join(row))
    print()

# Function to move a car
def move_car(b, car, direction, steps):
    # Find the car's position
    positions = [(r, c) for r in range(6) for c in range(6) if b[r][c] == car]
    if not positions:
        return False

    # Determine if the car is horizontal or vertical
    is_horizontal = positions[0][0] == positions[1][0]

    # Move the car
    if is_horizontal:
        # Horizontal move
        for r, c in positions:
            b[r][c] = '.'
        for r, c in positions:
            new_c = c + steps if direction == 'right' else c - steps
            b[r][new_c] = car
    else:
        # Vertical move
        for r, c in positions:
            b[r][c] = '.'
        for r, c in positions:
            new_r = r + steps if direction == 'down' else r - steps
            b[new_r][c] = car

    return True

# Sequence of moves to solve the puzzle
moves = [
    ('L', 'down', 1),  # Move L down
    ('K', 'down', 1),  # Move K down
    ('A', 'right', 2)  # Move A to the exit
]

# Apply the moves
for car, direction, steps in moves:
    move_car(board, car, direction, steps)

# Print the final board
print_board(board)