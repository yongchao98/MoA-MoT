# Define the initial state of the puzzle
initial_state = [
    ['+', '+', '+', '+', '+', '+'],
    ['+', '-', '-', '-', '-', '+'],
    ['+', '-', 'X', '-', '-', '+'],
    ['+', '@', '@', '-', '-', '+'],
    ['+', '%', '@', 'X', '$', '+'],
    ['+', '+', '+', '+', '+', '+']
]

# Define the target state of the puzzle
target_state = [
    ['+', '+', '+', '+', '+', '+'],
    ['+', '-', '-', '-', '-', '+'],
    ['+', '-', 'X', '-', '-', '+'],
    ['+', '-', '-', '-', '-', '+'],
    ['+', '-', '$', 'X', '$', '+'],
    ['+', '+', '+', '+', '+', '+']
]

# Define a new sequence of moves
moves = "UURRDLDDLUURR"

# Function to apply moves to the initial state
def apply_moves(state, moves):
    # Find the player's initial position
    for i, row in enumerate(state):
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
        if state[new_pos[0]][new_pos[1]] in ['-', 'X']:
            # Move the player
            state[i][j] = '-' if state[i][j] == '%' else 'X'
            state[new_pos[0]][new_pos[1]] = '%' if state[new_pos[0]][new_pos[1]] == 'X' else '*'
            player_pos = new_pos
        elif state[new_pos[0]][new_pos[1]] in ['@', '$']:
            # Check if the box can be moved
            box_new_pos = (new_pos[0] + (new_pos[0] - i), new_pos[1] + (new_pos[1] - j))
            if state[box_new_pos[0]][box_new_pos[1]] in ['-', 'X']:
                # Move the box
                state[new_pos[0]][new_pos[1]] = '-' if state[new_pos[0]][new_pos[1]] == '@' else 'X'
                state[box_new_pos[0]][box_new_pos[1]] = '@' if state[box_new_pos[0]][box_new_pos[1]] == '-' else '$'
                # Move the player
                state[i][j] = '-' if state[i][j] == '%' else 'X'
                state[new_pos[0]][new_pos[1]] = '%' if state[new_pos[0]][new_pos[1]] == 'X' else '*'
                player_pos = new_pos

    return state

# Apply the moves and check if the target state is reached
final_state = apply_moves(initial_state, moves)
is_solved = final_state == target_state

print("Is the puzzle solved?", is_solved)