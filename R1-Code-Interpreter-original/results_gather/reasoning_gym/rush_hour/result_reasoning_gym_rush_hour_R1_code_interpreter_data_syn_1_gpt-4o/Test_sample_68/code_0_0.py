# Initial board setup
board = [
    ['E', '.', 'F', 'B', 'B', 'J'],
    ['E', '.', 'F', '.', 'I', 'J'],
    ['E', 'A', 'A', 'G', 'I', 'K'],
    ['.', '.', '.', 'G', '.', 'K'],
    ['C', 'C', '.', 'H', '.', 'K'],
    ['D', 'D', 'D', 'H', '.', '.']
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
        for r, c in positions:
            b[r][c] = '.'
        for r, c in positions:
            if direction == 'right':
                b[r][c + steps] = car
            elif direction == 'left':
                b[r][c - steps] = car
    else:
        for r, c in positions:
            b[r][c] = '.'
        for r, c in positions:
            if direction == 'down':
                b[r + steps][c] = car
            elif direction == 'up':
                b[r - steps][c] = car

    return True

# Simulate the moves
print("Initial Board:")
print_board(board)

# Move K down
move_car(board, 'K', 'down', 1)
print("After moving K down:")
print_board(board)

# Move G down
move_car(board, 'G', 'down', 1)
print("After moving G down:")
print_board(board)

# Move I down
move_car(board, 'I', 'down', 1)
print("After moving I down:")
print_board(board)

# Move A to the right
move_car(board, 'A', 'right', 3)
print("After moving A to the right:")
print_board(board)

# Final moves
moves = "K+1 G+1 I+1 A+3"
print(moves)