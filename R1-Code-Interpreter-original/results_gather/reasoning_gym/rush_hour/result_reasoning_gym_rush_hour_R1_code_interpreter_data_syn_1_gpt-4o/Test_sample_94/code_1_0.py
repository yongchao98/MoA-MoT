# Initial board setup
board = [
    "..D...",
    "..DBBB",
    ".AAE.F",
    "...E.F",
    ".CCx.F",
    "......"
]

# Function to print the board
def print_board(board):
    for row in board:
        print(row)
    print()

# Function to move a car
def move_car(board, car, direction, steps):
    new_board = [list(row) for row in board]
    car_positions = [(i, j) for i, row in enumerate(board) for j, cell in enumerate(row) if cell == car]
    
    if not car_positions:
        return board  # Car not found, return the original board

    # Determine if the car is horizontal or vertical
    is_horizontal = car_positions[0][0] == car_positions[-1][0]

    # Clear the car's current position
    for i, j in car_positions:
        new_board[i][j] = '.'

    # Move the car
    for i, j in car_positions:
        if direction == 'right' and is_horizontal:
            new_board[i][j + steps] = car
        elif direction == 'left' and is_horizontal:
            new_board[i][j - steps] = car
        elif direction == 'down' and not is_horizontal:
            new_board[i + steps][j] = car
        elif direction == 'up' and not is_horizontal:
            new_board[i - steps][j] = car

    return [''.join(row) for row in new_board]

# Move sequence
moves = []

# Move F to the right by 1
board = move_car(board, 'F', 'right', 1)
moves.append('F+1')

# Move E down by 1
board = move_car(board, 'E', 'down', 1)
moves.append('E+1')

# Move AA to the right by 2
board = move_car(board, 'A', 'right', 2)
moves.append('A+2')

# Print the final board
print_board(board)

# Print the moves
print('Moves:', ' '.join(moves))