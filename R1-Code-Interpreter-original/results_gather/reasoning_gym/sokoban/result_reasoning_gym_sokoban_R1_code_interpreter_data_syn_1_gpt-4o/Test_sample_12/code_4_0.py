from heapq import heappop, heappush

def solve_sokoban():
    # Initial state of the puzzle
    initial_state = [
        "+ + + + + + + + + +",
        "+ $ + + + X - @ - +",
        "+ + + + + X @ @ - +",
        "+ X + + X - - X - +",
        "+ - + X @ - - X - +",
        "+ @ - - X @ - - - +",
        "+ - - - @ - @ - - +",
        "+ - - - - - - * - +",
        "+ + + + + + + + + +"
    ]

    # Convert the initial state to a more manageable format
    board = [list(row.replace(" ", "")) for row in initial_state]

    # Define the player's initial position
    player_pos = (7, 7)

    # Define the moves as (dy, dx) for U, D, L, R
    moves = {
        'U': (-1, 0),
        'D': (1, 0),
        'L': (0, -1),
        'R': (0, 1)
    }

    # A function to check if a move is valid
    def is_valid_move(board, player_pos, move):
        y, x = player_pos
        dy, dx = moves[move]
        new_y, new_x = y + dy, x + dx

        # Check if the new position is within bounds and not a wall
        if board[new_y][new_x] == '+':
            return False

        # Check if the player is pushing a box
        if board[new_y][new_x] in ('@', '$'):
            # Check if the box can be pushed
            box_new_y, box_new_x = new_y + dy, new_x + dx
            if board[box_new_y][box_new_x] in ('+', '@', '$'):
                return False

        return True

    # A function to apply a move
    def apply_move(board, player_pos, move):
        y, x = player_pos
        dy, dx = moves[move]
        new_y, new_x = y + dy, x + dx

        # Move the player
        if board[new_y][new_x] in ('-', 'X'):
            board[y][x] = '-' if board[y][x] == '*' else 'X'
            board[new_y][new_x] = '*' if board[new_y][new_x] == '-' else '%'
        elif board[new_y][new_x] in ('@', '$'):
            # Move the box
            box_new_y, box_new_x = new_y + dy, new_x + dx
            board[box_new_y][box_new_x] = '$' if board[box_new_y][box_new_x] == 'X' else '@'
            board[new_y][new_x] = '*' if board[new_y][new_x] == '@' else '%'
            board[y][x] = '-' if board[y][x] == '*' else 'X'

        return (new_y, new_x)

    # A function to check if the puzzle is solved
    def is_solved(board):
        for row in board:
            if '@' in row:
                return False
        return True

    # A heuristic function to estimate the cost to the goal
    def heuristic(board):
        # Simple heuristic: count the number of boxes not on goals
        return sum(row.count('@') for row in board)

    # A function to perform an A* search to find the solution
    def a_star(board, player_pos):
        open_set = []
        heappush(open_set, (0, board, player_pos, ""))
        visited = set()

        while open_set:
            _, current_board, current_pos, path = heappop(open_set)

            if is_solved(current_board):
                return path

            # Convert the board to a tuple of strings to use as a hashable state
            state = tuple(''.join(row) for row in current_board)
            if state in visited:
                continue
            visited.add(state)

            for move in moves:
                if is_valid_move(current_board, current_pos, move):
                    new_board = [row[:] for row in current_board]
                    new_player_pos = apply_move(new_board, current_pos, move)
                    cost = len(path) + 1 + heuristic(new_board)
                    heappush(open_set, (cost, new_board, new_player_pos, path + move))

        return None

    # Find the solution
    solution = a_star(board, player_pos)
    return solution

# Print the solution
print(solve_sokoban())