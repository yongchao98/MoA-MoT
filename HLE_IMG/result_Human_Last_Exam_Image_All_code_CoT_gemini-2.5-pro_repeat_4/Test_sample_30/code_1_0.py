import collections

def solve_knights_puzzle():
    """
    Solves the Knights Puzzle for the five given configurations.
    This function encapsulates the entire logic and prints the final result.
    """

    # --- Problem Setup ---

    # The 4x3 board is indexed as follows:
    # 0  1  2
    # 3  4  5
    # 6  7  8
    # 9 10 11
    MOVES = {
        0: [5, 7], 1: [6, 8], 2: [3, 7],
        3: [2, 8, 10], 4: [9, 11], 5: [0, 6, 10],
        6: [1, 5, 11], 7: [0, 2], 8: [1, 3, 9],
        9: [4, 8], 10: [3, 5], 11: [4, 6]
    }

    # Initial configurations from the image
    # Format: (frozenset_of_white_knight_positions, frozenset_of_black_knight_positions)
    CONFIGS = {
        'A': (frozenset([2, 5, 8, 11]), frozenset([0, 3, 6, 9])),
        'B': (frozenset([3, 9, 11]), frozenset([0, 6, 8])),
        'C': (frozenset([0, 7]), frozenset([1, 5])),
        'D': (frozenset([0, 7]), frozenset([4, 10])),
        'E': (frozenset([1, 2, 5]), frozenset([0, 3, 4]))
    }

    def is_solvable(initial_w, initial_b):
        """
        Determines if a configuration is solvable using Breadth-First Search.
        """
        # The goal is to swap the initial positions.
        goal_w = initial_b
        goal_b = initial_w

        # A state is (white_positions, black_positions, player_to_move).
        # We start with White's turn.
        start_state = (initial_w, initial_b, 'W')

        # queue for BFS, and a set to keep track of visited states to avoid cycles.
        queue = collections.deque([start_state])
        visited = {start_state}

        while queue:
            current_w, current_b, turn = queue.popleft()

            # Check if the current board configuration is the goal.
            if current_w == goal_w and current_b == goal_b:
                return True

            # Determine whose knights to move and whose are stationary.
            moving_knights = current_w if turn == 'W' else current_b
            other_knights = current_b if turn == 'W' else current_w
            all_occupied_squares = current_w | current_b
            next_turn = 'B' if turn == 'W' else 'W'

            # Generate all possible next states from the current state.
            for knight_pos in moving_knights:
                for move_dest in MOVES[knight_pos]:
                    # A move is valid only if the destination square is empty.
                    if move_dest not in all_occupied_squares:
                        new_moving_knights = (moving_knights - {knight_pos}) | {move_dest}

                        if turn == 'W':
                            new_w, new_b = new_moving_knights, other_knights
                        else:
                            new_w, new_b = other_knights, new_moving_knights

                        new_state = (new_w, new_b, next_turn)

                        if new_state not in visited:
                            visited.add(new_state)
                            queue.append(new_state)
        
        # If the queue becomes empty, all reachable states have been explored
        # and the goal was not found.
        return False

    # --- Main Execution ---
    solvable_configs = []
    print("Analyzing Knights Puzzle configurations...")
    print("-" * 30)
    for name, (w_pos, b_pos) in sorted(CONFIGS.items()):
        is_solvable_result = is_solvable(w_pos, b_pos)
        if is_solvable_result:
            solvable_configs.append(name)
        status = "Solvable" if is_solvable_result else "Unsolvable"
        print(f"Configuration {name}: {status}")
    
    print("-" * 30)
    print("Conclusion:")
    if not solvable_configs:
        print("None of the configurations are solvable.")
    else:
        # The final answer format is specific, so we print the solvable list.
        print(f"The solvable configurations are: {', '.join(solvable_configs)}")

# Run the solver
solve_knights_puzzle()