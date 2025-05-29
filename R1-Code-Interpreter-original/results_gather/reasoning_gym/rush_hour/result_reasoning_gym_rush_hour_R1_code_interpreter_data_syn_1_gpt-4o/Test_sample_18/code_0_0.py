# Initial board setup
board = [
    ['B', 'B', 'H', 'I', 'C', 'C'],
    ['.', 'G', 'H', 'I', 'J', '.'],
    ['.', 'G', 'A', 'A', 'J', 'K'],
    ['D', 'D', '.', '.', 'x', 'K'],
    ['E', 'E', '.', '.', '.', 'K'],
    ['.', '.', 'x', 'F', 'F', '.']
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
    is_horizontal = positions[0][0] == positions[1][0]

    # Move the car
    if is_horizontal:
        if direction == '+':
            # Move right
            for r, c in reversed(positions):
                b[r][c + spaces] = car
                b[r][c] = '.'
        elif direction == '-':
            # Move left
            for r, c in positions:
                b[r][c - spaces] = car
                b[r][c] = '.'
    else:
        if direction == '+':
            # Move down
            for r, c in reversed(positions):
                b[r + spaces][c] = car
                b[r][c] = '.'
        elif direction == '-':
            # Move up
            for r, c in positions:
                b[r - spaces][c] = car
                b[r][c] = '.'

    return True

# Simulate the moves
print("Initial Board:")
print_board(board)

# Move J down by 1
move_car(board, 'J', '+', 1)
print("After moving J+1:")
print_board(board)

# Move K down by 1
move_car(board, 'K', '+', 1)
print("After moving K+1:")
print_board(board)

# Move AA to the right by 2
move_car(board, 'A', '+', 2)
print("After moving A+2:")
print_board(board)