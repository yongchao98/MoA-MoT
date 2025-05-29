# Initial board setup
board = [
    ['B', 'B', 'F', '.', '.', 'H'],
    ['.', 'E', 'F', '.', '.', 'H'],
    ['.', 'E', 'F', 'A', 'A', 'I'],
    ['.', 'C', 'C', 'G', '.', 'I'],
    ['x', '.', '.', 'G', 'D', 'D'],
    ['.', '.', '.', 'x', '.', '.']
]

# Function to print the board
def print_board(board):
    for row in board:
        print(''.join(row))
    print()

# Function to move a car
def move_car(board, car, direction, spaces):
    if direction == 'right':
        for i in range(len(board)):
            for j in range(len(board[i]) - 1, -1, -1):
                if board[i][j] == car:
                    board[i][j] = '.'
                    board[i][j + spaces] = car
    elif direction == 'down':
        for i in range(len(board) - 1, -1, -1):
            for j in range(len(board[i])):
                if board[i][j] == car:
                    board[i][j] = '.'
                    board[i + spaces][j] = car

# Simulate the moves
print("Initial Board:")
print_board(board)

# Move H down by 2
move_car(board, 'H', 'down', 2)
print("After moving H down by 2:")
print_board(board)

# Move I down by 1
move_car(board, 'I', 'down', 1)
print("After moving I down by 1:")
print_board(board)

# Move AA to the right by 2
move_car(board, 'A', 'right', 2)
print("After moving AA to the right by 2:")
print_board(board)