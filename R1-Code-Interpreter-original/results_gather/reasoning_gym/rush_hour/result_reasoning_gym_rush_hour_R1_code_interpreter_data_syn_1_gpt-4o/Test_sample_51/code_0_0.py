# Initial board setup
board = [
    ['B', 'B', 'B', 'K', 'C', 'C'],
    ['D', 'D', 'J', 'K', '.', 'L'],
    ['I', '.', 'J', 'A', 'A', 'L'],
    ['I', 'E', 'E', '.', 'x', 'M'],
    ['F', 'F', '.', '.', '.', 'M'],
    ['G', 'G', '.', 'x', 'H', 'H']
]

# Function to print the board
def print_board(b):
    for row in b:
        print(''.join(row))
    print()

# Move L down by 1
board[1][5], board[2][5] = '.', 'L'

# Move AA right by 1
board[2][3], board[2][4], board[2][5] = '.', 'A', 'A'

# Print the final board
print_board(board)