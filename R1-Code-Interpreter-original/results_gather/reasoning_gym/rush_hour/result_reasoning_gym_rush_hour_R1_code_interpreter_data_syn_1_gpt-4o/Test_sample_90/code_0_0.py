# Initial board setup
board = [
    ['B', 'B', 'B', 'K', 'C', 'C'],
    ['D', 'D', 'J', 'K', '.', 'L'],
    ['I', '.', 'J', 'A', 'A', 'L'],
    ['I', 'E', 'E', '.', 'x', 'M'],
    ['F', 'F', '.', '.', '.', 'M'],
    ['G', 'G', 'H', 'H', '.', 'x']
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
        row = positions[0][0]
        if direction == 'right':
            for _ in range(steps):
                if positions[-1][1] + 1 < 6 and b[row][positions[-1][1] + 1] == '.':
                    for r, c in reversed(positions):
                        b[r][c + 1] = car
                        b[r][c] = '.'
                    positions = [(r, c + 1) for r, c in positions]
                else:
                    return False
        elif direction == 'left':
            for _ in range(steps):
                if positions[0][1] - 1 >= 0 and b[row][positions[0][1] - 1] == '.':
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
                if positions[-1][0] + 1 < 6 and b[positions[-1][0] + 1][col] == '.':
                    for r, c in reversed(positions):
                        b[r + 1][c] = car
                        b[r][c] = '.'
                    positions = [(r + 1, c) for r, c in positions]
                else:
                    return False
        elif direction == 'up':
            for _ in range(steps):
                if positions[0][0] - 1 >= 0 and b[positions[0][0] - 1][col] == '.':
                    for r, c in positions:
                        b[r - 1][c] = car
                        b[r][c] = '.'
                    positions = [(r - 1, c) for r, c in positions]
                else:
                    return False
    return True

# Sequence of moves to solve the puzzle
moves = [
    ('L', 'down', 2),  # Move L down 2 spaces
    ('A', 'right', 1)  # Move A right 1 space to the exit
]

# Apply the moves
for car, direction, steps in moves:
    move_car(board, car, direction, steps)

# Print the final board
print_board(board)