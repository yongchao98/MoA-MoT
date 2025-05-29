# Initial board setup
board = [
    ['x', '.', 'B', 'B', 'J', '.'],
    ['.', '.', 'H', '.', 'J', 'K'],
    ['G', '.', 'H', 'A', 'A', 'K'],
    ['G', 'C', 'C', 'I', '.', 'K'],
    ['G', 'D', 'D', 'I', 'E', 'E'],
    ['x', '.', 'F', 'F', '.', '.']
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
            b[r][c + steps] = car
    else:
        for r, c in positions:
            b[r][c] = '.'
        for r, c in positions:
            b[r + steps][c] = car

    return True

# Simulate the moves
print("Initial board:")
print_board(board)

# Move H left by 1
move_car(board, 'H', 'left', -1)
print("After moving H left by 1:")
print_board(board)

# Move J up by 1
move_car(board, 'J', 'up', -1)
print("After moving J up by 1:")
print_board(board)

# Move K down by 1
move_car(board, 'K', 'down', 1)
print("After moving K down by 1:")
print_board(board)

# Move AA right by 1 to exit
move_car(board, 'A', 'right', 1)
print("After moving AA right by 1 to exit:")
print_board(board)

# Final moves
moves = "H-1 J-1 K+1 A+1"
print(moves)