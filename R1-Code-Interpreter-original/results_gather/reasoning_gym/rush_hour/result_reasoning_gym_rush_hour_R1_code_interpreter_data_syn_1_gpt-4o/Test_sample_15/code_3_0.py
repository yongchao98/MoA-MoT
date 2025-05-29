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
    car_positions = []
    for i, row in enumerate(board):
        for j, cell in enumerate(row):
            if cell == car:
                car_positions.append((i, j))

    if not car_positions:
        return  # Car not found

    # Determine the orientation of the car
    if len(car_positions) > 1 and car_positions[0][0] == car_positions[1][0]:
        orientation = 'horizontal'
    else:
        orientation = 'vertical'

    # Move the car
    if orientation == 'horizontal':
        car_row = car_positions[0][0]
        if direction == 'right':
            for _ in range(steps):
                if car_positions[-1][1] + 1 < len(board[0]) and board[car_row][car_positions[-1][1] + 1] == '.':
                    # Move the car to the right
                    board[car_row] = board[car_row][:car_positions[0][1]] + '.' * len(car_positions) + board[car_row][car_positions[-1][1] + 1:]
                    car_positions = [(car_row, j + 1) for _, j in car_positions]
                    for _, j in car_positions:
                        board[car_row] = board[car_row][:j] + car + board[car_row][j + 1:]
                else:
                    break
    else:
        car_col = car_positions[0][1]
        if direction == 'down':
            for _ in range(steps):
                if car_positions[-1][0] + 1 < len(board) and board[car_positions[-1][0] + 1][car_col] == '.':
                    # Move the car down
                    for i, j in car_positions:
                        board[i] = board[i][:j] + '.' + board[i][j + 1:]
                    car_positions = [(i + 1, car_col) for i, _ in car_positions]
                    for i, j in car_positions:
                        board[i] = board[i][:j] + car + board[i][j + 1:]
                else:
                    break

# Simulate the moves
print("Initial board:")
print_board(board)

# Move L down
move_car(board, 'L', 'down', 1)
print("After moving L down:")
print_board(board)

# Move J down
move_car(board, 'J', 'down', 1)
print("After moving J down:")
print_board(board)

# Move K down
move_car(board, 'K', 'down', 1)
print("After moving K down:")
print_board(board)

# Move A to the right
move_car(board, 'A', 'right', 2)
print("After moving A to the right:")
print_board(board)