# Define the initial board state
board = [
    ['.', 'B', 'B', 'I', 'C', 'C'],
    ['x', 'D', 'D', 'I', '.', 'J'],
    ['x', '.', 'H', 'A', 'A', 'J'],
    ['G', '.', 'H', '.', '.', 'K'],
    ['G', 'E', 'E', '.', '.', 'K'],
    ['.', '.', '.', 'F', 'F', 'F']
]

# Define the moves to solve the puzzle
moves = [
    ('J', 1),  # Move J down by 1
    ('H', -1), # Move H up by 1
    ('A', 2)   # Move A right by 2
]

# Function to apply moves to the board
def apply_moves(board, moves):
    for car, move in moves:
        if car == 'J':
            # Move J down by 1
            board[1][5], board[2][5] = '.', 'J'
        elif car == 'H':
            # Move H up by 1
            board[2][2], board[1][2] = '.', 'H'
        elif car == 'A':
            # Move A right by 2
            board[2][3], board[2][4], board[2][5] = '.', '.', 'A'
    return board

# Apply the moves and print the final board state
final_board = apply_moves(board, moves)
for row in final_board:
    print(''.join(row))