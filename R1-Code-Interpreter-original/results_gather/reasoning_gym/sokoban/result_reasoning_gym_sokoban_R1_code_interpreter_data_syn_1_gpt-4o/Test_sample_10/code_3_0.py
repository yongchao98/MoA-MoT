def solve_sokoban():
    # Initial state of the puzzle
    puzzle = [
        ['+', '+', '+', '+', '+', '+'],
        ['+', 'X', 'X', '@', '-', '+'],
        ['+', '-', 'X', '-', '-', '+'],
        ['+', '@', '-', '@', '-', '+'],
        ['+', '%', '@', '-', '-', '+'],
        ['+', '+', '+', '+', '+', '+']
    ]

    # Find the player's initial position
    player_pos = (4, 1)

    # Moves: L, R, U, D
    moves = "UURRDLUURRDDLL"

    # Directions corresponding to moves
    directions = {
        'L': (0, -1),
        'R': (0, 1),
        'U': (-1, 0),
        'D': (1, 0)
    }

    # Apply the moves to the puzzle
    for move in moves:
        dx, dy = directions[move]
        new_x, new_y = player_pos[0] + dx, player_pos[1] + dy

        # Check if the new position is a wall
        if puzzle[new_x][new_y] == '+':
            continue

        # Check if the new position is a box
        if puzzle[new_x][new_y] in ('@', '$'):
            # Calculate the position behind the box
            box_new_x, box_new_y = new_x + dx, new_y + dy

            # Check if the box can be moved
            if puzzle[box_new_x][box_new_y] in ('-', 'X'):
                # Move the box
                if puzzle[box_new_x][box_new_y] == 'X':
                    puzzle[box_new_x][box_new_y] = '$'
                else:
                    puzzle[box_new_x][box_new_y] = '@'

                # Update the current box position
                if puzzle[new_x][new_y] == '$':
                    puzzle[new_x][new_y] = 'X'
                else:
                    puzzle[new_x][new_y] = '-'

        # Move the player
        if puzzle[new_x][new_y] == 'X':
            puzzle[new_x][new_y] = '%'
        else:
            puzzle[new_x][new_y] = '@'

        # Update the player's previous position
        if puzzle[player_pos[0]][player_pos[1]] == '%':
            puzzle[player_pos[0]][player_pos[1]] = 'X'
        else:
            puzzle[player_pos[0]][player_pos[1]] = '-'

        # Update the player's position
        player_pos = (new_x, new_y)

    # Return the solution moves
    return moves

# Print the solution
print(solve_sokoban())