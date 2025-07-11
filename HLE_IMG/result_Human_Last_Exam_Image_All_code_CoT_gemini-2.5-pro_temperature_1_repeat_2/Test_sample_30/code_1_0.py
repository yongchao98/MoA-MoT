import collections

def is_solvable(initial_board):
    """
    Checks if a given Knights Puzzle configuration is solvable using Breadth-First Search.

    A state is defined by a tuple: (board_configuration, current_turn).
    The board configuration is a tuple of strings, making it hashable.
    The turn is 'W' for White or 'B' for Black.
    """
    
    # 1. Define the goal state by swapping W and B on the initial board.
    goal_list = []
    for row in initial_board:
        new_row = ""
        for char in row:
            if char == 'W':
                new_row += 'B'
            elif char == 'B':
                new_row += 'W'
            else:
                new_row += '.'
        goal_list.append(new_row)
    goal_board = tuple(goal_list)
    
    # 2. Initialize the BFS queue and visited set.
    # The queue stores states to visit: (board, turn).
    # White moves first.
    q = collections.deque([(initial_board, 'W')])
    # The visited set prevents cycles and redundant explorations.
    visited = {(initial_board, 'W')}
    
    # All 8 possible L-shaped knight moves.
    moves = [(-2, -1), (-2, 1), (-1, -2), (-1, 2), 
             (1, -2), (1, 2), (2, -1), (2, 1)]

    # 3. Start the BFS loop.
    while q:
        current_board, current_turn = q.popleft()

        # Check if the current board configuration is the goal.
        if current_board == goal_board:
            return True

        # Determine which pieces to move and what the next turn will be.
        piece_to_move = 'W' if current_turn == 'W' else 'B'
        next_turn = 'B' if current_turn == 'W' else 'W'

        # 4. Generate next possible states.
        # Iterate through the board to find all knights of the current player.
        for r in range(4): # 4 rows
            for c in range(3): # 3 columns
                if current_board[r][c] == piece_to_move:
                    # For each knight, try all 8 possible moves.
                    for dr, dc in moves:
                        nr, nc = r + dr, c + dc

                        # Check if the move is valid: within bounds and to an empty square.
                        if 0 <= nr < 4 and 0 <= nc < 3 and current_board[nr][nc] == '.':
                            # Create the new board configuration after the move.
                            new_board_list = [list(row) for row in current_board]
                            new_board_list[r][c] = '.'
                            new_board_list[nr][nc] = piece_to_move
                            new_board = tuple("".join(row) for row in new_board_list)
                            
                            # The new state includes the new board and the next player's turn.
                            new_state = (new_board, next_turn)

                            # If this new state has not been visited, add it to the queue.
                            if new_state not in visited:
                                visited.add(new_state)
                                q.append(new_state)

    # 5. If the queue is empty and the goal was not found, the puzzle is unsolvable.
    return False

def main():
    """
    Defines the five configurations from the image and solves them.
    """
    # Define the initial configurations. 'B' = Black, 'W' = White, '.' = Empty.
    configurations = {
        'A': (
            "B.W",
            "B.W",
            "B.W",
            "B.W"
        ),
        'B': (
            ".B.",
            ".W.",
            "B.B",
            "W.W"
        ),
        'C': (
            "W.B",
            "..B",
            ".W.",
            "..."
        ),
        'D': (
            "W..",
            ".B.",
            ".W.",
            ".B."
        ),
        'E': (
            "BWW",
            "BBW",
            "...",
            "..."
        )
    }

    solvable_configs = []
    for name, config in configurations.items():
        if is_solvable(config):
            solvable_configs.append(name)
    
    solvable_configs.sort()
    print("The following configurations are solvable:")
    print(", ".join(solvable_configs))

if __name__ == "__main__":
    main()