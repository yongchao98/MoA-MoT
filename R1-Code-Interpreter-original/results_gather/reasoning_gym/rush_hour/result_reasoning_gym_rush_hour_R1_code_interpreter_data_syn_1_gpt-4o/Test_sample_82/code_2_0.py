# Initial board setup
board = [
    ['F', '.', '.', 'H', 'B', 'B'],
    ['F', 'C', 'C', 'H', '.', 'I'],
    ['A', 'A', 'G', '.', '.', 'I'],
    ['.', '.', 'G', 'D', 'D', 'I'],
    ['.', '.', 'x', 'x', '.', '.'],
    ['.', '.', '.', 'E', 'E', 'E']
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

    # Move the car
    if is_horizontal:
        # Horizontal move
        for r, c in positions:
            b[r][c] = '.'
        for r, c in positions:
            if direction == '+':
                b[r][c + steps] = car
            else:
                b[r][c - steps] = car
    else:
        # Vertical move
        for r, c in positions:
            b[r][c] = '.'
        for r, c in positions:
            if direction == '+':
                b[r + steps][c] = car
            else:
                b[r - steps][c] = car

    return True

# Simulate the moves
print("Initial Board:")
print_board(board)

# Move I up
move_car(board, 'I', '-', 1)
print("After moving I up:")
print_board(board)

# Move G up
move_car(board, 'G', '-', 1)
print("After moving G up:")
print_board(board)

# Move AA to the exit
move_car(board, 'A', '+', 3)
print("After moving AA to the exit:")
print_board(board)