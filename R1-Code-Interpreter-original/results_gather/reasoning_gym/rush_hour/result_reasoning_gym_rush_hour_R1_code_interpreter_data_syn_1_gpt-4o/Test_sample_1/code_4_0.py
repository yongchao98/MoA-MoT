# Initial board setup
board = [
    ['.', 'F', '.', 'I', 'B', 'B'],
    ['.', 'F', '.', 'I', '.', 'J'],
    ['.', 'G', 'A', 'A', '.', 'J'],
    ['.', 'G', 'C', 'C', 'C', 'J'],
    ['.', 'G', 'H', 'D', 'D', 'D'],
    ['.', '.', 'H', '.', 'E', 'E']
]

# Function to print the board
def print_board(b):
    for row in b:
        print(''.join(row))
    print()

# Function to move a car
def move_car(b, car, direction, steps):
    # Find the car's position
    positions = [(i, j) for i in range(6) for j in range(6) if b[i][j] == car]
    if not positions:
        return False

    # Determine if the car is vertical or horizontal
    vertical = all(p[1] == positions[0][1] for p in positions)

    # Move the car
    if vertical:
        # Vertical car
        if direction == 'up':
            for _ in range(steps):
                if positions[0][0] == 0 or b[positions[0][0] - 1][positions[0][1]] != '.':
                    return False
                for p in positions:
                    b[p[0] - 1][p[1]] = car
                    b[p[0]][p[1]] = '.'
                positions = [(p[0] - 1, p[1]) for p in positions]
        elif direction == 'down':
            for _ in range(steps):
                if positions[-1][0] == 5 or b[positions[-1][0] + 1][positions[-1][1]] != '.':
                    return False
                for p in reversed(positions):
                    b[p[0] + 1][p[1]] = car
                    b[p[0]][p[1]] = '.'
                positions = [(p[0] + 1, p[1]) for p in positions]
    else:
        # Horizontal car
        if direction == 'left':
            for _ in range(steps):
                if positions[0][1] == 0 or b[positions[0][0]][positions[0][1] - 1] != '.':
                    return False
                for p in positions:
                    b[p[0]][p[1] - 1] = car
                    b[p[0]][p[1]] = '.'
                positions = [(p[0], p[1] - 1) for p in positions]
        elif direction == 'right':
            for _ in range(steps):
                if positions[-1][1] == 5 or b[positions[-1][0]][positions[-1][1] + 1] != '.':
                    return False
                for p in reversed(positions):
                    b[p[0]][p[1] + 1] = car
                    b[p[0]][p[1]] = '.'
                positions = [(p[0], p[1] + 1) for p in positions]

    return True

# Sequence of moves to solve the puzzle
moves = [
    ('C', 'left', 1),
    ('J', 'down', 1),
    ('A', 'right', 2)
]

# Apply the moves
for car, direction, steps in moves:
    if not move_car(board, car, direction, steps):
        print(f"Failed to move {car} {direction} by {steps} steps.")
    print_board(board)

# Final board state
print_board(board)