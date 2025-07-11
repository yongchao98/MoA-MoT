import collections

def solve_knights_puzzle():
    """
    Solves the Knights Puzzle for configurations A-E on a 4x3 board.
    This function models the puzzle as a state-space search and uses BFS
    to determine if a solution exists for each configuration.
    """

    # Board dimensions and adjacency list for knight moves
    rows, cols = 4, 3
    adj = {
        0: [5, 7], 1: [6, 8], 2: [3, 7],
        3: [2, 8, 10], 4: [9, 11], 5: [0, 6, 10],
        6: [1, 5, 11], 7: [0, 2], 8: [1, 3],
        9: [4, 8], 10: [3, 5], 11: [4, 6]
    }

    # Define initial configurations from the image
    # A: B at {0,3,6,9}, W at {2,5,8,11}
    # B: B at {1,5,6}, W at {4,9,11}
    # C: B at {2,5}, W at {0,7}
    # D: B at {4,10}, W at {0,6}
    # E: B at {0,1,3,4}, W at {2,5}
    configs = {
        'A': ({2, 5, 8, 11}, {0, 3, 6, 9}),
        'B': ({4, 9, 11}, {1, 5, 6}),
        'C': ({0, 7}, {2, 5}),
        'D': ({0, 6}, {4, 10}),
        'E': ({2, 5}, {0, 1, 3, 4})
    }

    solvable_configs = []

    for name, (white_pos, black_pos) in configs.items():
        print(f"Analyzing configuration {name}...")
        
        if len(white_pos) != len(black_pos):
            print(f"Result: Configuration {name} is unsolvable (unequal number of knights).")
            continue

        # Define initial and target states as tuples
        board = ['E'] * (rows * cols)
        for pos in white_pos:
            board[pos] = 'W'
        for pos in black_pos:
            board[pos] = 'B'
        initial_state = tuple(board)

        target_board = ['E'] * (rows * cols)
        for pos in black_pos:
            target_board[pos] = 'W'
        for pos in white_pos:
            target_board[pos] = 'B'
        target_state = tuple(target_board)

        # BFS setup: queue stores (state, player_to_move)
        queue = collections.deque([(initial_state, 'W')])
        visited = {initial_state}
        
        is_solvable = False
        
        # Limit search to prevent potential infinite loops on very complex/large boards
        max_states_to_check = 500000 
        count = 0

        while queue:
            count += 1
            if count > max_states_to_check:
                print(f"Result: Configuration {name} is unsolvable (search limit exceeded).")
                is_solvable = False # Assumed unsolvable
                break

            current_state, player_to_move = queue.popleft()

            if current_state == target_state:
                is_solvable = True
                break

            # Generate next states by moving the current player's pieces
            for pos_idx, piece in enumerate(current_state):
                if piece == player_to_move:
                    for move_to in adj[pos_idx]:
                        if current_state[move_to] == 'E':
                            new_state_list = list(current_state)
                            new_state_list[pos_idx] = 'E'
                            new_state_list[move_to] = player_to_move
                            new_state = tuple(new_state_list)

                            if new_state not in visited:
                                visited.add(new_state)
                                next_player = 'B' if player_to_move == 'W' else 'W'
                                queue.append((new_state, next_player))
        
        if is_solvable:
            print(f"Result: Configuration {name} is solvable.")
            solvable_configs.append(name)
        elif count <= max_states_to_check:
             print(f"Result: Configuration {name} is unsolvable (entire state space searched).")


    print("\n--- Summary ---")
    if solvable_configs:
        solvable_str = ", ".join(sorted(solvable_configs))
        print(f"The solvable configurations are: {solvable_str}")
    else:
        print("None of the configurations are solvable.")

# Run the solver
solve_knights_puzzle()