# Initial board setup
board = [
    ['.', '.', '.', 'G', '.', '.'],
    ['.', '.', '.', 'G', 'B', 'B'],
    ['.', '.', 'A', 'A', 'H', 'I'],
    ['x', '.', '.', '.', 'H', 'I'],
    ['F', 'C', 'C', 'D', 'D', 'I'],
    ['F', '.', '.', 'E', 'E', 'E']
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

# Initial board
print("Initial board:")
print_board(board)

# Move 'H' down by 1
move_car(board, 'H', 'down', 1)
print("After moving H down by 1:")
print_board(board)

# Move 'I' down by 1
move_car(board, 'I', 'down', 1)
print("After moving I down by 1:")
print_board(board)

# Move 'A' right by 2
move_car(board, 'A', 'right', 2)
print("After moving A right by 2:")
print_board(board)

# Check if 'A' is at the exit
a_at_exit = board[2][4] == 'A' and board[2][5] == 'A'
print("A at exit:", a_at_exit)