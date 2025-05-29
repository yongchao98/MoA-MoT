# Define the initial board state
board = [
    ['B', 'B', 'H', '.', 'C', 'C'],
    ['G', '.', 'H', '.', 'J', 'K'],
    ['G', 'A', 'A', '.', 'J', 'K'],
    ['D', 'D', '.', 'I', 'x', 'L'],
    ['E', 'E', '.', 'I', '.', 'L'],
    ['F', 'F', 'F', '.', '.', 'x']
]

# Define the moves
moves = [
    ('H', 1),  # Move H down
    ('J', 1),  # Move J down
    ('K', 1),  # Move K down
    ('A', 3)   # Move A to the right to the exit
]

# Apply the moves to the board
for car, move in moves:
    if car == 'A':
        # Move the red car to the right
        board[2][1] = '.'
        board[2][2] = '.'
        board[2][3] = 'A'
        board[2][4] = 'A'
        board[2][5] = 'A'
    elif car == 'H':
        # Move H down
        board[0][2] = '.'
        board[1][2] = 'H'
        board[2][2] = 'H'
    elif car == 'J':
        # Move J down
        board[1][4] = '.'
        board[2][4] = 'J'
        board[3][4] = 'J'
    elif car == 'K':
        # Move K down
        board[1][5] = '.'
        board[2][5] = 'K'
        board[3][5] = 'K'

# Print the final board state
for row in board:
    print(''.join(row))