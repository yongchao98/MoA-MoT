# Initial board setup
board = [
    ['G', '.', 'I', 'B', 'B', 'x'],
    ['G', '.', 'I', 'J', 'C', 'C'],
    ['H', 'A', 'A', 'J', '.', 'K'],
    ['H', 'D', 'D', 'D', '.', 'K'],
    ['E', 'E', 'F', 'F', '.', 'K'],
    ['.', '.', '.', '.', '.', '.']
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
                    board[i][j + spaces] = car
                    board[i][j] = '.'
    elif direction == 'down':
        for i in range(len(board) - 1, -1, -1):
            for j in range(len(board[i])):
                if board[i][j] == car:
                    board[i + spaces][j] = car
                    board[i][j] = '.'

# Initial board
print("Initial board:")
print_board(board)

# Move J down by 1
move_car(board, 'J', 'down', 1)
print("After moving J down by 1:")
print_board(board)

# Move K down by 1
move_car(board, 'K', 'down', 1)
print("After moving K down by 1:")
print_board(board)

# Move AA to the right by 3
move_car(board, 'A', 'right', 3)
print("After moving AA to the right by 3:")
print_board(board)