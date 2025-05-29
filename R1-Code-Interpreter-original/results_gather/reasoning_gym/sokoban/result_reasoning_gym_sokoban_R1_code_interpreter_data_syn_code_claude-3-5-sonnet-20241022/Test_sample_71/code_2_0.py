def solve_sokoban(board):
    def get_player_pos():
        for i in range(len(board)):
            for j in range(len(board[i])):
                if board[i][j] in ['*', '%']:
                    return (i, j)
        return None

    def is_valid(r, c):
        return 0 <= r < len(board) and 0 <= c < len(board[0]) and board[r][c] != '+'

    def get_board_state():
        state = []
        for i in range(len(board)):
            row = []
            for j in range(len(board[i])):
                row.append(board[i][j])
            state.append(row)
        return state

    def set_cell(r, c, val):
        board[r][c] = val

    def try_move(path, depth=0):
        if depth > 50:  # Limit search depth
            return None
            
        # Check if solved
        boxes_left = 0
        goals_unfilled = 0
        for row in board:
            for cell in row:
                if cell == '@':
                    boxes_left += 1
                if cell == 'X':
                    goals_unfilled += 1
        if boxes_left == 0 and goals_unfilled == 0:
            return path

        player = get_player_pos()
        if not player:
            return None

        pr, pc = player
        # Try each direction: Left, Up, Right, Down
        directions = [('L', 0, -1), ('U', -1, 0), ('R', 0, 1), ('D', 1, 0)]
        
        for move, dr, dc in directions:
            new_r, new_c = pr + dr, pc + dc
            
            if not is_valid(new_r, new_c):
                continue

            # Save current state
            old_state = get_board_state()
            moved = False

            # Moving to empty space or goal
            if board[new_r][new_c] in ['-', 'X']:
                if board[pr][pc] == '*':
                    set_cell(pr, pc, 'X')
                else:
                    set_cell(pr, pc, '-')
                
                if board[new_r][new_c] == 'X':
                    set_cell(new_r, new_c, '%')
                else:
                    set_cell(new_r, new_c, '*')
                moved = True

            # Pushing a box
            elif board[new_r][new_c] in ['@', '$']:
                next_r, next_c = new_r + dr, new_c + dc
                if is_valid(next_r, next_c) and board[next_r][next_c] in ['-', 'X']:
                    # Move box
                    if board[next_r][next_c] == 'X':
                        set_cell(next_r, next_c, '$')
                    else:
                        set_cell(next_r, next_c, '@')
                    
                    # Move player
                    if board[pr][pc] == '*':
                        set_cell(pr, pc, 'X')
                    else:
                        set_cell(pr, pc, '-')
                    
                    if board[new_r][new_c] == '$':
                        set_cell(new_r, new_c, '%')
                    else:
                        set_cell(new_r, new_c, '*')
                    moved = True

            if moved:
                result = try_move(path + move, depth + 1)
                if result:
                    return result
                    
                # Restore state
                for i in range(len(board)):
                    for j in range(len(board[i])):
                        board[i][j] = old_state[i][j]

        return None

    return try_move("")

# Initialize the board
board = [
    ['+', '+', '+', '+', '+', '+', '+'],
    ['+', '-', '-', '-', '-', '-', '+'],
    ['+', '-', '-', '@', '-', '-', '+'],
    ['+', '-', '-', '-', '-', 'X', '+'],
    ['+', 'X', '@', '-', 'X', '-', '+'],
    ['+', '-', '@', '*', '@', '-', '+'],
    ['+', '-', '-', 'X', '@', 'X', '+'],
    ['+', '+', '+', '+', '+', '+', '+']
]

solution = solve_sokoban(board)
print(solution)