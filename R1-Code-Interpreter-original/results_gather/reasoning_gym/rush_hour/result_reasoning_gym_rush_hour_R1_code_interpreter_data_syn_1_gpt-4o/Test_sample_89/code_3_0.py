# Initial board setup
board = [
    ['I', 'B', 'B', 'C', 'C', '.'],
    ['I', 'D', 'D', 'E', 'E', 'L'],
    ['I', 'A', 'A', '.', 'K', 'L'],
    ['F', 'F', 'F', '.', 'K', '.'],
    ['.', 'J', 'G', 'G', 'H', 'H'],
    ['.', 'J', '.', '.', '.', '.']
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
print("Initial board:")
print_board(board)

# Move L down by 1
move_car(board, 'L', 'down', 1)
print("After moving L down by 1:")
print_board(board)

# Move K down by 1
move_car(board, 'K', 'down', 1)
print("After moving K down by 1:")
print_board(board)

# Move AA right by 4
move_car(board, 'A', 'right', 4)
print("After moving AA right by 4:")
print_board(board)