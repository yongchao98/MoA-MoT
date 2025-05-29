# Initial board setup
board = [
    ['.', '.', 'x', 'B', 'B', 'L'],
    ['.', '.', 'I', 'C', 'C', 'L'],
    ['A', 'A', 'I', 'J', '.', '.'],
    ['H', 'D', 'D', 'J', 'E', 'E'],
    ['H', 'F', 'F', 'J', 'K', 'x'],
    ['G', 'G', '.', '.', 'K', '.']
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
        return False

    # Determine if the car is horizontal or vertical
    is_horizontal = positions[0][0] == positions[1][0]

    # Move the car
    if is_horizontal:
        if direction == 'right':
            for r, c in reversed(positions):
                if c + steps >= len(b[0]) or b[r][c + steps] not in ('.', car):
                    return False
                b[r][c] = '.'
                b[r][c + steps] = car
        elif direction == 'left':
            for r, c in positions:
                if c - steps < 0 or b[r][c - steps] not in ('.', car):
                    return False
                b[r][c] = '.'
                b[r][c - steps] = car
    else:
        if direction == 'down':
            for r, c in reversed(positions):
                if r + steps >= len(b) or b[r + steps][c] not in ('.', car):
                    return False
                b[r][c] = '.'
                b[r + steps][c] = car
        elif direction == 'up':
            for r, c in positions:
                if r - steps < 0 or b[r - steps][c] not in ('.', car):
                    return False
                b[r][c] = '.'
                b[r - steps][c] = car

    return True

# Simulate the moves
moves = [
    ('J', 'down', 1),
    ('K', 'right', 1),
    ('A', 'right', 3)
]

# Apply the moves
for car, direction, steps in moves:
    if not move_car(board, car, direction, steps):
        print(f"Failed to move {car} {direction} by {steps} steps.")
        break

# Print the final board
print_board(board)