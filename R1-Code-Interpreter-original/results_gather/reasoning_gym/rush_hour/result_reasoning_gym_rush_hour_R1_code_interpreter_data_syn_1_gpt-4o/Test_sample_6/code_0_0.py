# Initial board setup
board = [
    ['B', 'B', 'H', '.', 'K', '.'],
    ['.', 'G', 'H', '.', 'K', 'L'],
    ['.', 'G', 'A', 'A', 'K', 'L'],
    ['C', 'C', 'I', '.', '.', 'x'],
    ['F', '.', 'I', 'J', 'D', 'D'],
    ['F', 'E', 'E', 'J', '.', 'x']
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
                if c + steps < len(b[0]) and b[r][c + steps] == '.':
                    b[r][c + steps] = car
                    b[r][c] = '.'
                else:
                    return False
        elif direction == 'left':
            for r, c in positions:
                if c - steps >= 0 and b[r][c - steps] == '.':
                    b[r][c - steps] = car
                    b[r][c] = '.'
                else:
                    return False
    else:
        if direction == 'down':
            for r, c in reversed(positions):
                if r + steps < len(b) and b[r + steps][c] == '.':
                    b[r + steps][c] = car
                    b[r][c] = '.'
                else:
                    return False
        elif direction == 'up':
            for r, c in positions:
                if r - steps >= 0 and b[r - steps][c] == '.':
                    b[r - steps][c] = car
                    b[r][c] = '.'
                else:
                    return False
    return True

# Perform the moves
moves = [
    ('K', 'down', 1),
    ('L', 'down', 1),
    ('A', 'right', 2)
]

# Apply the moves
for car, direction, steps in moves:
    if not move_car(board, car, direction, steps):
        print(f"Move {car} {direction} {steps} failed.")
        break

# Print the final board
print_board(board)