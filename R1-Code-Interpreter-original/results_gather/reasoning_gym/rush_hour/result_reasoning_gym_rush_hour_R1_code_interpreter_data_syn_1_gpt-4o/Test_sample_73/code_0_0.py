# Initial board setup
board = [
    ['H', 'I', 'B', 'B', 'C', 'C'],
    ['H', 'I', 'J', '.', '.', '.'],
    ['A', 'A', 'J', '.', '.', 'L'],
    ['D', 'D', 'J', 'K', '.', 'L'],
    ['x', 'E', 'E', 'K', '.', 'L'],
    ['.', 'F', 'F', 'G', 'G', '.']
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
        row = positions[0][0]
        if direction == 'right':
            for _ in range(steps):
                if b[row][positions[-1][1] + 1] == '.':
                    for r, c in reversed(positions):
                        b[r][c + 1] = car
                        b[r][c] = '.'
                    positions = [(r, c + 1) for r, c in positions]
                else:
                    return False
        elif direction == 'left':
            for _ in range(steps):
                if b[row][positions[0][1] - 1] == '.':
                    for r, c in positions:
                        b[r][c - 1] = car
                        b[r][c] = '.'
                    positions = [(r, c - 1) for r, c in positions]
                else:
                    return False
    else:
        col = positions[0][1]
        if direction == 'down':
            for _ in range(steps):
                if b[positions[-1][0] + 1][col] == '.':
                    for r, c in reversed(positions):
                        b[r + 1][c] = car
                        b[r][c] = '.'
                    positions = [(r + 1, c) for r, c in positions]
                else:
                    return False
        elif direction == 'up':
            for _ in range(steps):
                if b[positions[0][0] - 1][col] == '.':
                    for r, c in positions:
                        b[r - 1][c] = car
                        b[r][c] = '.'
                    positions = [(r - 1, c) for r, c in positions]
                else:
                    return False
    return True

# Simulate the moves
moves = [
    ('J', 'down', 1),
    ('K', 'down', 1),
    ('L', 'down', 1),
    ('A', 'right', 3)
]

# Apply the moves
for car, direction, steps in moves:
    if not move_car(board, car, direction, steps):
        print(f"Failed to move {car} {direction} by {steps} steps.")
        break

# Print the final board
print_board(board)