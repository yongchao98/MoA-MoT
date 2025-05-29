# Initial board setup
board = [
    "..xBBB",
    "GHI.CC",
    "GHIAAK",
    "DDDJ.K",
    "...JEE",
    ".FFF.x"
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
                elif direction == 'down':
                    new_board[i + steps][j] = car
    return [''.join(row) for row in new_board]

# Initial board
print("Initial Board:")
print_board(board)

# Move K down by 1
board = move_car(board, 'K', 'down', 1)
print("After moving K down by 1:")
print_board(board)

# Move AA right by 1 to the exit
board = move_car(board, 'A', 'right', 1)
print("After moving AA right by 1 to the exit:")
print_board(board)

# Final moves
moves = "K+1 A+1"
print(moves)