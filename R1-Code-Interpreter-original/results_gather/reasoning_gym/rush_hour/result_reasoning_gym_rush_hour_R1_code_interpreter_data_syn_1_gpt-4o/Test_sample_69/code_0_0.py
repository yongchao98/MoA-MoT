# Initial board setup
board = [
    ['G', 'B', 'B', 'B', 'I', '.'],
    ['G', 'C', 'C', '.', 'I', '.'],
    ['A', 'A', 'H', '.', 'I', 'J'],
    ['.', 'x', 'H', 'D', 'D', 'J'],
    ['E', 'E', 'E', '.', '.', 'J'],
    ['.', '.', 'F', 'F', 'F', '.']
]

# Moves to solve the puzzle
moves = []

# Move H down by 1
board[2][2], board[3][2] = '.', 'H'
moves.append('H+1')

# Move AA to the right by 4
board[2][0], board[2][1], board[2][2], board[2][3], board[2][4] = '.', '.', 'A', 'A', '.'
moves.append('A+4')

# Print the moves
print(' '.join(moves))