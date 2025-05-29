# Initial board setup
board = [
    ['F', 'B', 'B', 'C', 'C', 'x'],
    ['F', '.', 'G', '.', 'D', 'D'],
    ['A', 'A', 'G', '.', 'H', 'I'],
    ['.', '.', 'G', '.', 'H', 'I'],
    ['.', '.', '.', '.', '.', 'J'],
    ['E', 'E', 'E', '.', '.', 'J']
]

# Moves to solve the puzzle
moves = [
    ('G', 1),  # Move G down
    ('H', 1),  # Move H down
    ('I', 1),  # Move I down
    ('J', 1),  # Move J down
    ('D', 1),  # Move D right
    ('C', 1),  # Move C right
    ('A', 4)   # Move A right to the exit
]

# Function to apply moves
def apply_moves(board, moves):
    for car, spaces in moves:
        if car in ['A', 'B', 'C', 'D', 'E', 'F']:  # Horizontal cars
            for i in range(len(board)):
                for j in range(len(board[i])):
                    if board[i][j] == car:
                        # Move the car horizontally
                        board[i][j] = '.'
                        board[i][j + spaces] = car
                        break
        else:  # Vertical cars
            for j in range(len(board[0])):
                for i in range(len(board)):
                    if board[i][j] == car:
                        # Move the car vertically
                        board[i][j] = '.'
                        board[i + spaces][j] = car
                        break
    return board

# Apply the moves
final_board = apply_moves(board, moves)

# Print the final board state
for row in final_board:
    print(''.join(row))