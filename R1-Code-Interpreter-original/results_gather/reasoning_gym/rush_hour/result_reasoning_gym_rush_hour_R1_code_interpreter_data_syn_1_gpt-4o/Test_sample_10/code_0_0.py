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
    for i, row in enumerate(board):
        for j, cell in enumerate(row):
            if cell == car:
                new_board[i][j] = '.'
                if direction == 'right':
                    new_board[i][j + steps] = car
                elif direction == 'left':
                    new_board[i][j - steps] = car
                elif direction == 'down':
                    new_board[i + steps][j] = car
                elif direction == 'up':
                    new_board[i - steps][j] = car
    return ["".join(row) for row in new_board]

# Move K down by 1
board = move_car(board, 'K', 'down', 1)
print_board(board)

# Move A right by 3
board = move_car(board, 'A', 'right', 3)
print_board(board)

# The moves are: K+1 A+3
print("<<<K+1 A+3>>>")