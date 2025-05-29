# Define the initial state of the board
board = [
    ['+', '+', '+', '+', '+', '+'],
    ['+', '-', '$', '+', '$', '+'],
    ['+', '-', '-', '-', '+', '+'],
    ['+', '-', '@', '-', '+', '+'],
    ['+', '-', '-', 'X', 'X', '+'],
    ['+', '-', '@', 'X', '@', '+'],
    ['+', '-', '-', '-', '*', '+'],
    ['+', '+', '+', '+', '+', '+']
]

# Define the player's initial position
player_pos = (6, 4)

# Define the sequence of moves to solve the puzzle
# This is a hypothetical solution and may need adjustments
moves = "UULLDRRURRDDLLURR"

# Function to apply moves to the board
def apply_moves(board, player_pos, moves):
    for move in moves:
        x, y = player_pos
        if move == 'U':
            new_pos = (x - 1, y)
        elif move == 'D':
            new_pos = (x + 1, y)
        elif move == 'L':
            new_pos = (x, y - 1)
        elif move == 'R':
            new_pos = (x, y + 1)
        
        # Check if the new position is valid
        if board[new_pos[0]][new_pos[1]] in ['-', 'X']:
            # Move the player
            board[x][y] = '-'
            board[new_pos[0]][new_pos[1]] = '*'
            player_pos = new_pos
        elif board[new_pos[0]][new_pos[1]] in ['@', '$']:
            # Check if the box can be pushed
            box_new_pos = (new_pos[0] + (new_pos[0] - x), new_pos[1] + (new_pos[1] - y))
            if board[box_new_pos[0]][box_new_pos[1]] in ['-', 'X']:
                # Move the box
                board[new_pos[0]][new_pos[1]] = '*'
                board[box_new_pos[0]][box_new_pos[1]] = '@' if board[box_new_pos[0]][box_new_pos[1]] == '-' else '$'
                board[x][y] = '-'
                player_pos = new_pos
    return board

# Apply the moves to the board
final_board = apply_moves(board, player_pos, moves)

# Print the final board state
for row in final_board:
    print(' '.join(row))