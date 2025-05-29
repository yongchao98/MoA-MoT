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
    new_positions = []
    for i, j in car_positions:
        if direction == 'right':
            new_j = j + steps
            if new_j >= len(board[0]) or any(new_board[i][j + k] != '.' for k in range(1, steps + 1)):
                return board  # Invalid move
            new_positions.append((i, new_j))
        elif direction == 'left':
            new_j = j - steps
            if new_j < 0 or any(new_board[i][j - k] != '.' for k in range(1, steps + 1)):
                return board  # Invalid move
            new_positions.append((i, new_j))
        elif direction == 'down':
            new_i = i + steps
            if new_i >= len(board) or any(new_board[i + k][j] != '.' for k in range(1, steps + 1)):
                return board  # Invalid move
            new_positions.append((new_i, j))
        elif direction == 'up':
            new_i = i - steps
            if new_i < 0 or any(new_board[i - k][j] != '.' for k in range(1, steps + 1)):
                return board  # Invalid move
            new_positions.append((new_i, j))

    # Clear the car's current positions
    for i, j in car_positions:
        new_board[i][j] = '.'

    # Move the car to the new positions
    for i, j in new_positions:
        new_board[i][j] = car

    return ["".join(row) for row in new_board]

# Move K down by 1
board = move_car(board, 'K', 'down', 1)
print_board(board)

# Move A right by 2
board = move_car(board, 'A', 'right', 2)
print_board(board)

# The moves are: K+1 A+2
print("<<<K+1 A+2>>>")