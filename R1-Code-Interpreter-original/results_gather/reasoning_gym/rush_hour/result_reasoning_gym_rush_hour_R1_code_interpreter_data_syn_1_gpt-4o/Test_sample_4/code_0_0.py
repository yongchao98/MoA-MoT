# Initial board setup
board = [
    ['B', 'B', 'G', 'I', '.', 'K'],
    ['E', 'F', 'G', 'I', 'x', 'K'],
    ['E', 'F', 'A', 'A', 'J', 'L'],
    ['C', 'C', 'H', '.', 'J', 'L'],
    ['.', '.', 'H', 'D', 'D', '.'],
    ['.', '.', '.', '.', '.', '.']
]

# Function to print the board
def print_board(b):
    for row in b:
        print(''.join(row))
    print()

# Move J down by 1
board[2][4], board[3][4] = '.', 'J'

# Move AA right by 2
board[2][2], board[2][3], board[2][4], board[2][5] = '.', '.', 'A', 'A'

# Print the final board
print_board(board)