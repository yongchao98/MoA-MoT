# Initial board setup
board = [
    ['.', 'B', 'B', 'C', 'C', 'J'],
    ['x', '.', '.', '.', '.', 'J'],
    ['G', 'A', 'A', 'I', '.', 'K'],
    ['G', 'D', 'D', 'I', '.', 'K'],
    ['.', '.', 'H', 'I', 'E', 'E'],
    ['F', 'F', 'H', '.', '.', '.']
]

# Function to print the board
def print_board(b):
    for row in b:
        print(''.join(row))
    print()

# Function to move a car
def move_car(b, car, direction, steps):
    # Find the car's position
    positions = [(r, c) for r in range(len(b)) for c in range(len(b[0])) if b[r][c] == car]
    if not positions:
        return False  # Car not found

    # Determine if the car is horizontal or vertical
    is_horizontal = positions[0][0] == positions[1][0]

    # Move the car
    if is_horizontal:
        if direction == 'right':
            for r, c in reversed(positions):
                if c + steps >= len(b[0]) or b[r][c + steps] not in ('.', car):
                    return False  # Move is invalid
                b[r][c + steps] = car
                b[r][c] = '.'
        elif direction == 'left':
            for r, c in positions:
                if c - steps < 0 or b[r][c - steps] not in ('.', car):
                    return False  # Move is invalid
                b[r][c - steps] = car
                b[r][c] = '.'
    else:
        if direction == 'down':
            for r, c in reversed(positions):
                if r + steps >= len(b) or b[r + steps][c] not in ('.', car):
                    return False  # Move is invalid
                b[r + steps][c] = car
                b[r][c] = '.'
        elif direction == 'up':
            for r, c in positions:
                if r - steps < 0 or b[r - steps][c] not in ('.', car):
                    return False  # Move is invalid
                b[r - steps][c] = car
                b[r][c] = '.'

    return True

# Perform the moves
moves = [
    ('I', 'down', 1),
    ('K', 'down', 1),
    ('A', 'right', 3)
]

# Apply the moves
for car, direction, steps in moves:
    if not move_car(board, car, direction, steps):
        print(f"Move {car} {direction} {steps} failed.")
    print_board(board)

# Final board state
print_board(board)