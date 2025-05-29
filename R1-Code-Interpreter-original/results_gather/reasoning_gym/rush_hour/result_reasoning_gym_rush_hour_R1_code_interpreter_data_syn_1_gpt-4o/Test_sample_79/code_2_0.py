# Initial board setup
board = [
    ['.', 'B', 'B', 'C', 'C', '.'],
    ['x', 'x', 'D', 'D', 'K', '.'],
    ['G', 'A', 'A', 'J', 'K', '.'],
    ['G', 'E', 'E', 'J', '.', '.'],
    ['.', 'H', 'I', 'J', '.', '.'],
    ['.', 'H', 'I', 'F', 'F', 'F']
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

    # Calculate new positions
    new_positions = []
    for r, c in positions:
        if is_horizontal:
            new_c = c + steps if direction in ['right', 'left'] else c
            new_positions.append((r, new_c))
        else:
            new_r = r + steps if direction in ['down', 'up'] else r
            new_positions.append((new_r, c))

    # Check if new positions are valid
    for r, c in new_positions:
        if r < 0 or r >= len(b) or c < 0 or c >= len(b[0]) or (b[r][c] != '.' and b[r][c] != car):
            return False  # Invalid move

    # Move the car
    for r, c in positions:
        b[r][c] = '.'
    for r, c in new_positions:
        b[r][c] = car

    return True

# Simulate the moves
print("Initial board:")
print_board(board)

# Move K down
move_car(board, 'K', 'down', 1)
print("After moving K down:")
print_board(board)

# Move E down
move_car(board, 'E', 'down', 1)
print("After moving E down:")
print_board(board)

# Move A to the right
move_car(board, 'A', 'right', 2)
print("After moving A to the right:")
print_board(board)

# Output the moves
print("Moves: <<<K+1 E+1 A+2>>>")