# Initial board setup
board = [
    ['x', 'B', 'B', 'B', '.', 'I'],
    ['.', 'C', 'C', '.', '.', 'I'],
    ['.', '.', 'G', 'A', 'A', 'I'],
    ['.', '.', 'G', 'H', 'D', 'D'],
    ['.', '.', 'G', 'H', 'E', 'E'],
    ['.', 'F', 'F', 'x', '.', '.']
]

# Function to print the board
def print_board(b):
    for row in b:
        print(''.join(row))
    print()

# Function to move a car
def move_car(b, car, direction, spaces):
    # Find the car's position
    positions = [(r, c) for r in range(len(b)) for c in range(len(b[0])) if b[r][c] == car]
    if not positions:
        return False

    # Determine the orientation of the car
    vertical = all(p[1] == positions[0][1] for p in positions)
    horizontal = all(p[0] == positions[0][0] for p in positions)

    # Calculate new positions
    new_positions = []
    if vertical:
        for r, c in positions:
            new_r = r + spaces if direction == 'down' else r - spaces
            new_positions.append((new_r, c))
    elif horizontal:
        for r, c in positions:
            new_c = c + spaces if direction == 'right' else c - spaces
            new_positions.append((r, new_c))

    # Check if new positions are valid
    for r, c in new_positions:
        if r < 0 or r >= len(b) or c < 0 or c >= len(b[0]) or (b[r][c] != '.' and b[r][c] != car):
            return False

    # Move the car
    for r, c in positions:
        b[r][c] = '.'
    for r, c in new_positions:
        b[r][c] = car

    return True

# Simulate the moves
print("Initial Board:")
print_board(board)

# Move I up
move_car(board, 'I', 'up', 1)
print("After moving I up:")
print_board(board)

# Move D down
move_car(board, 'D', 'down', 1)
print("After moving D down:")
print_board(board)

# Move A to the right
move_car(board, 'A', 'right', 2)
print("After moving A to the right:")
print_board(board)