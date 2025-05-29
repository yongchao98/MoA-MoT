# Initial board setup
board = [
    ".xBBCC",
    "..x.K.",
    "G.AAK.",
    "G.IJDD",
    "H.IJ..",
    "HEEFFF"
]

# Function to print the board
def print_board(board):
    for row in board:
        print(row)
    print()

# Function to move a car
def move_car(board, car, direction, steps):
    new_board = [list(row) for row in board]
    car_positions = []

    # Find all positions of the car
    for i, row in enumerate(board):
        for j, cell in enumerate(row):
            if cell == car:
                car_positions.append((i, j))
    
    # Determine the new positions and check boundaries
    if direction == 'right':
        # Check if the path is clear for the entire car to move right
        for i, j in car_positions:
            if j + steps >= len(board[0]) or new_board[i][j + steps] != '.':
                return board  # Invalid move
        # Move the car to the new positions
        for i, j in reversed(car_positions):
            new_board[i][j] = '.'
            new_board[i][j + steps] = car
    elif direction == 'left':
        # Check if the path is clear for the entire car to move left
        for i, j in car_positions:
            if j - steps < 0 or new_board[i][j - steps] != '.':
                return board  # Invalid move
        # Move the car to the new positions
        for i, j in car_positions:
            new_board[i][j] = '.'
            new_board[i][j - steps] = car
    elif direction == 'down':
        # Check if the path is clear for the entire car to move down
        for i, j in car_positions:
            if i + steps >= len(board) or new_board[i + steps][j] != '.':
                return board  # Invalid move
        # Move the car to the new positions
        for i, j in reversed(car_positions):
            new_board[i][j] = '.'
            new_board[i + steps][j] = car
    elif direction == 'up':
        # Check if the path is clear for the entire car to move up
        for i, j in car_positions:
            if i - steps < 0 or new_board[i - steps][j] != '.':
                return board  # Invalid move
        # Move the car to the new positions
        for i, j in car_positions:
            new_board[i][j] = '.'
            new_board[i - steps][j] = car

    return ["".join(row) for row in new_board]

# Move K down by 1
board = move_car(board, 'K', 'down', 1)
print_board(board)

# Move A right by 2
board = move_car(board, 'A', 'right', 2)
print_board(board)

# The moves are: K+1 A+2
print("<<<K+1 A+2>>>")