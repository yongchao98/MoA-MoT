# Initial board setup
board = [
    ['B', 'B', 'B', 'H', '.', '.'],
    ['F', '.', '.', 'H', '.', '.'],
    ['F', 'A', 'A', 'I', '.', 'J'],
    ['C', 'C', 'G', 'I', '.', 'J'],
    ['.', '.', 'G', 'D', 'D', 'J'],
    ['x', 'x', 'E', 'E', '.', '.']
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
        if direction == '+':  # Move right
            for _ in range(steps):
                if positions[-1][1] + 1 < len(b[0]) and b[row][positions[-1][1] + 1] == '.':
                    for r, c in reversed(positions):
                        b[r][c + 1] = car
                        b[r][c] = '.'
                    positions = [(r, c + 1) for r, c in positions]
                else:
                    return False
        elif direction == '-':  # Move left
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
        if direction == '+':  # Move down
            for _ in range(steps):
                if positions[-1][0] + 1 < len(b) and b[positions[-1][0] + 1][col] == '.':
                    for r, c in reversed(positions):
                        b[r + 1][c] = car
                        b[r][c] = '.'
                    positions = [(r + 1, c) for r, c in positions]
                else:
                    return False
        elif direction == '-':  # Move up
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
    ('J', '+', 1),  # Move J down
    ('I', '+', 1),  # Move I right
    ('A', '+', 4)   # Move AA right to exit
]

# Apply the moves
for car, direction, steps in moves:
    move_car(board, car, direction, steps)

# Print the final board
print_board(board)