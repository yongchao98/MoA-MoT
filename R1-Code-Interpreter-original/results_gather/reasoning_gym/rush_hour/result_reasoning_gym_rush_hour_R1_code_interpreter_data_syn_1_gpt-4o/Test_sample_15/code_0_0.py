# Initial board setup
board = [
    ".IBBBx",
    "GI.JCC",
    "GAAJK.",
    "G.DDKL",
    "H.xEEL",
    "HFFF.."
]

# Function to print the board
def print_board(board):
    for row in board:
        print(row)
    print()

# Function to move a car
def move_car(board, car, direction, steps):
    # Find the car on the board
    for i, row in enumerate(board):
        if car in row:
            car_row = i
            car_col = row.index(car)
            break

    # Determine the orientation of the car
    if board[car_row][car_col + 1] == car:
        orientation = 'horizontal'
    else:
        orientation = 'vertical'

    # Move the car
    if orientation == 'horizontal':
        if direction == 'right':
            for _ in range(steps):
                if board[car_row][car_col + 2] == '.':
                    board[car_row] = board[car_row][:car_col] + '.' + car + car + board[car_row][car_col + 3:]
                    car_col += 1
                else:
                    break
    else:
        if direction == 'down':
            for _ in range(steps):
                if board[car_row + 2][car_col] == '.':
                    board[car_row] = board[car_row][:car_col] + '.' + board[car_row][car_col + 1:]
                    board[car_row + 2] = board[car_row + 2][:car_col] + car + board[car_row + 2][car_col + 1:]
                    car_row += 1
                else:
                    break

# Simulate the moves
print("Initial board:")
print_board(board)

# Move J down
move_car(board, 'J', 'down', 1)
print("After moving J down:")
print_board(board)

# Move K down
move_car(board, 'K', 'down', 1)
print("After moving K down:")
print_board(board)

# Move L down
move_car(board, 'L', 'down', 1)
print("After moving L down:")
print_board(board)

# Move A to the right
move_car(board, 'A', 'right', 2)
print("After moving A to the right:")
print_board(board)