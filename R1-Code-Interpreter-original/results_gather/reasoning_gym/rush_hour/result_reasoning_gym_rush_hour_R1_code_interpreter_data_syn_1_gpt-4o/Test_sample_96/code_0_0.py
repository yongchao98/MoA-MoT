# Initial board setup
board = [
    ['B', 'B', 'C', 'C', '.', 'x'],
    ['D', 'D', 'J', 'E', 'E', 'M'],
    ['.', '.', 'J', 'A', 'A', 'M'],
    ['x', '.', 'K', 'F', 'F', 'M'],
    ['.', '.', 'K', 'L', 'G', 'G'],
    ['.', 'H', 'H', 'L', 'I', 'I']
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
    for r, c in positions:
        b[r][c] = '.'

    if is_horizontal:
        if direction == 'right':
            for r, c in positions:
                b[r][c + steps] = car
        elif direction == 'left':
            for r, c in positions:
                b[r][c - steps] = car
    else:
        if direction == 'down':
            for r, c in positions:
                b[r + steps][c] = car
        elif direction == 'up':
            for r, c in positions:
                b[r - steps][c] = car

    return True

# Sequence of moves
moves = [
    ('M', 'down', 1),  # Move M down
    ('A', 'right', 1)  # Move A to the exit
]

# Apply the moves
for car, direction, steps in moves:
    move_car(board, car, direction, steps)
    print_board(board)

# Final board state
print_board(board)