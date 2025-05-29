# Initial board setup
board = [
    ['B', 'B', 'B', 'F', '.', '.'],
    ['C', 'C', 'C', 'F', '.', '.'],
    ['.', 'A', 'A', 'F', '.', '.'],
    ['.', '.', 'E', '.', '.', '.'],
    ['.', '.', 'E', 'D', 'D', 'D'],
    ['.', '.', '.', '.', '.', '.']
]

# Function to print the board
def print_board(board):
    for row in board:
        print(''.join(row))
    print()

# Move car F down by 1
board[0][3], board[1][3], board[2][3], board[3][3] = '.', '.', 'F', 'F'

# Move car C right by 1
board[1][0], board[1][1], board[1][2], board[1][3] = '.', 'C', 'C', 'C'

# Move car F up by 1
board[2][3], board[3][3], board[4][3] = 'F', '.', 'F'

# Move car A right by 3
board[2][1], board[2][2], board[2][3], board[2][4], board[2][5] = '.', '.', '.', 'A', 'A'

# Print the final board
print_board(board)

# Output the moves
print("<<<F+1 C+1 F-1 A+3>>>")