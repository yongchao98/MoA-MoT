# Define the initial state of the puzzle
puzzle = [
    ['+', '+', '+', '+', '+', '+'],
    ['+', '-', '-', '@', 'X', '+'],
    ['+', '%', '@', '-', '-', '+'],
    ['+', '-', '-', '$', '-', '+'],
    ['+', '-', '-', '-', '-', '+'],
    ['+', '+', 'X', '@', '-', '+'],
    ['+', '+', '+', '+', '+', '+']
]

# Define the new moves
moves = "RRDLLDRR"

# Function to apply moves to the puzzle
def apply_moves(puzzle, moves):
    # Find the player's initial position
    for i, row in enumerate(puzzle):
        for j, cell in enumerate(row):
            if cell == '%' or cell == '*':
                player_pos = (i, j)
                break

    # Apply each move
    for move in moves:
        i, j = player_pos
        if move == 'U':
            new_pos = (i - 1, j)
        elif move == 'D':
            new_pos = (i + 1, j)
        elif move == 'L':
            new_pos = (i, j - 1)
        elif move == 'R':
            new_pos = (i, j + 1)
        
        # Check if the move is valid
        if puzzle[new_pos[0]][new_pos[1]] in ['-', 'X']:
            # Move the player
            puzzle[i][j] = '-' if puzzle[i][j] == '%' else 'X'
            puzzle[new_pos[0]][new_pos[1]] = '%' if puzzle[new_pos[0]][new_pos[1]] == '-' else '*'
            player_pos = new_pos
        elif puzzle[new_pos[0]][new_pos[1]] in ['@', '$']:
            # Check if the box can be moved
            box_new_pos = (new_pos[0] + (new_pos[0] - i), new_pos[1] + (new_pos[1] - j))
            if puzzle[box_new_pos[0]][box_new_pos[1]] in ['-', 'X']:
                # Move the box
                puzzle[new_pos[0]][new_pos[1]] = '%' if puzzle[new_pos[0]][new_pos[1]] == '@' else '*'
                puzzle[box_new_pos[0]][box_new_pos[1]] = '@' if puzzle[box_new_pos[0]][box_new_pos[1]] == '-' else '$'
                puzzle[i][j] = '-' if puzzle[i][j] == '%' else 'X'
                player_pos = new_pos

    return puzzle

# Apply the moves and print the final state
final_state = apply_moves(puzzle, moves)
for row in final_state:
    print(' '.join(row))