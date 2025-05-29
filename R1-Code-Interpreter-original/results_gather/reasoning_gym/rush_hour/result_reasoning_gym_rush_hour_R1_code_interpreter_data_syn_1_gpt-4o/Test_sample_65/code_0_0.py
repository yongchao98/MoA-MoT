# Initial board setup
board = [
    ['E', '.', '.', 'H', '.', '.'],
    ['E', '.', '.', 'H', 'x', 'x'],
    ['E', 'A', 'A', 'I', 'J', '.'],
    ['.', 'F', 'G', 'I', 'J', 'K'],
    ['.', 'F', 'G', 'B', 'B', 'K'],
    ['C', 'C', 'D', 'D', 'D', '.']
]

# Function to print the board
def print_board(b):
    for row in b:
        print(''.join(row))
    print()

# Move 'I' and 'J' down
board[2][3] = '.'
board[2][4] = '.'
board[3][3] = 'I'
board[3][4] = 'J'
board[4][3] = 'I'
board[4][4] = 'J'

# Move 'AA' to the right
board[2][2] = '.'
board[2][3] = 'A'
board[2][4] = 'A'

# Print the final board
print_board(board)