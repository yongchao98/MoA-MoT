# Initial board setup
board = [
    ['B', 'B', 'x', 'C', 'C', 'K'],
    ['G', '.', '.', 'I', '.', 'K'],
    ['G', 'A', 'A', 'I', '.', 'L'],
    ['.', '.', 'H', 'D', 'D', 'L'],
    ['.', '.', 'H', '.', 'J', '.'],
    ['E', 'E', 'F', 'F', 'J', 'x']
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
    for r, c in positions:
        b[r][c] = '.'

    if is_horizontal:
        if direction == '+':
            for r, c in positions:
                b[r][c + steps] = car
        else:
            for r, c in positions:
                b[r][c - steps] = car
    else:
        if direction == '+':
            for r, c in positions:
                b[r + steps][c] = car
        else:
            for r, c in positions:
                b[r - steps][c] = car

    return True

# Simulate the moves
print("Initial Board:")
print_board(board)

# Move I down by 1
move_car(board, 'I', '+', 1)
print("After moving I down by 1:")
print_board(board)

# Move L right by 1
move_car(board, 'L', '+', 1)
print("After moving L right by 1:")
print_board(board)

# Move A right by 3 to the exit
move_car(board, 'A', '+', 3)
print("After moving A right by 3 to the exit:")
print_board(board)

# Output the moves
print("<<<I+1 L+1 A+3>>>")