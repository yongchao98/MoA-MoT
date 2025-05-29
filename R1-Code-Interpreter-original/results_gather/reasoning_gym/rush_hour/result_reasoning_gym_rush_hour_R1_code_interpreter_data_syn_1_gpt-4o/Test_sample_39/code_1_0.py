# Initial board setup
board = [
    ['.', '.', 'x', 'B', 'B', 'K'],
    ['C', 'C', 'G', 'H', 'J', 'K'],
    ['A', 'A', 'G', 'H', 'J', '.'],
    ['.', '.', 'G', 'I', '.', '.'],
    ['F', 'D', 'D', 'I', '.', '.'],
    ['F', 'E', 'E', '.', '.', 'x']
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
        # Horizontal move
        for r, c in positions:
            b[r][c] = '.'
        for r, c in positions:
            new_c = c + steps if direction == 'right' else c - steps
            if 0 <= new_c < len(b[0]) and b[r][new_c] == '.':
                b[r][new_c] = car
            else:
                return False
    else:
        # Vertical move
        for r, c in positions:
            b[r][c] = '.'
        for r, c in positions:
            new_r = r + steps if direction == 'down' else r - steps
            if 0 <= new_r < len(b) and b[new_r][c] == '.':
                b[new_r][c] = car
            else:
                return False
    return True

# Sequence of moves to solve the puzzle
moves = [
    ('B', 'right', 1),
    ('K', 'right', 1),
    ('C', 'up', 1),
    ('G', 'up', 1),
    ('H', 'up', 1),
    ('J', 'up', 1),
    ('I', 'up', 1),
    ('D', 'left', 1),
    ('E', 'left', 1),
    ('A', 'right', 3)
]

# Apply the moves
for car, direction, steps in moves:
    move_car(board, car, direction, steps)

# Print the final board
print_board(board)