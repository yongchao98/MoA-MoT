# Initial board setup
board = [
    ['G', 'B', 'B', 'I', 'J', 'L'],
    ['G', '.', 'H', 'I', 'J', 'L'],
    ['A', 'A', 'H', '.', 'K', 'L'],
    ['C', 'C', 'D', 'D', 'K', '.'],
    ['.', '.', 'x', '.', '.', '.'],
    ['.', 'E', 'E', 'F', 'F', '.']
]

# Function to print the board
def print_board(b):
    for row in b:
        print(''.join(row))
    print()

# Move the car 'H' down by 1 space
board[1][2], board[2][2], board[3][2] = '.', '.', 'H'

# Move the car 'K' down by 1 space
board[2][4], board[3][4], board[4][4] = '.', '.', 'K'

# Move the red car 'AA' right by 3 spaces to reach the exit
board[2][0], board[2][1], board[2][2], board[2][3], board[2][4], board[2][5] = '.', '.', '.', '.', 'A', 'A'

# Print the final board
print_board(board)