# Initial board setup
board = [
    ['.', 'F', '.', 'H', 'B', 'B'],
    ['.', 'F', 'x', 'H', '.', 'J'],
    ['A', 'A', 'G', '.', '.', 'J'],
    ['E', '.', 'G', 'C', 'C', 'J'],
    ['E', '.', 'G', 'I', '.', '.'],
    ['E', 'D', 'D', 'I', '.', '.']
]

# Function to print the board
def print_board(b):
    for row in b:
        print(''.join(row))
    print()

# Move G down by 1
board[2][2], board[3][2], board[4][2] = '.', 'G', 'G'

# Move A to the right by 3
board[2][0], board[2][1], board[2][2] = '.', '.', 'A'
board[2][3], board[2][4] = 'A', '.'

# Print the final board
print_board(board)

# Output the moves
print("<<<G+1 A+3>>>")