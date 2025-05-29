# Define the initial board state
board = [
    ['.', 'H', 'B', 'B', 'C', 'C'],
    ['.', 'H', '.', '.', 'K', '.'],
    ['.', 'H', 'A', 'A', 'K', 'L'],
    ['D', 'D', 'I', '.', 'K', 'L'],
    ['G', '.', 'I', 'J', 'E', 'E'],
    ['G', 'F', 'F', 'J', '.', 'x']
]

# Define the moves
moves = [
    ('K', 1),  # Move K down by 1
    ('A', 2)   # Move A right by 2
]

# Apply the moves to the board
for car, move in moves:
    if car == 'K':
        # Move K down by 1
        board[1][4], board[2][4], board[3][4], board[4][4] = '.', 'K', 'K', 'K'
    elif car == 'A':
        # Move A right by 2
        board[2][2], board[2][3], board[2][4], board[2][5] = '.', '.', 'A', 'A'

# Print the final board state
for row in board:
    print(''.join(row))