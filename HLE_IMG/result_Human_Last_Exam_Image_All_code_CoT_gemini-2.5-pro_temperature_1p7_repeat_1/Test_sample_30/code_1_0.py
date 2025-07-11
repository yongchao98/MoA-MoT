import collections

def solve_knights_puzzle():
    """
    Solves the Knights Puzzle for five configurations on a 4x3 board.
    Determines which configurations are solvable.
    """
    # Board properties
    WIDTH, HEIGHT = 3, 4
    BOARD_SIZE = WIDTH * HEIGHT

    # Pre-calculate knight moves for each square (0-11)
    MOVES = {}
    for r_start in range(HEIGHT):
        for c_start in range(WIDTH):
            p_start = r_start * WIDTH + c_start
            MOVES[p_start] = []
            # All 8 possible knight move deltas
            for dr, dc in [(1, 2), (1, -2), (-1, 2), (-1, -2), (2, 1), (2, -1), (-2, 1), (-2, -1)]:
                r_end, c_end = r_start + dr, c_start + dc
                if 0 <= r_end < HEIGHT and 0 <= c_end < WIDTH:
                    p_end = r_end * WIDTH + c_end
                    MOVES[p_start].append(p_end)

    # Define initial configurations from the image
    # Board indexed 0-11, top-to-bottom, left-to-right.
    #   0  1  2
    #   3  4  5
    #   6  7  8
    #   9 10 11
    configs = {
        'A': {'black': [0, 3, 6, 9], 'white': [2, 5, 8, 11]},
        'B': {'black': [1, 4, 6, 8], 'white': [3, 9, 11]}, # 4 black, 3 white
        'C': {'black': [2, 5], 'white': [0, 7]},
        'D': {'black': [1, 5, 10], 'white': [0, 6, 7]},
        'E': {'black': [0, 3, 4], 'white': [1, 2, 5]}
    }
    
    solvable_configs = []

    print("Analyzing Knights Puzzle Configurations...\n")

    for name, pos in configs.items():
        initial_white_pos = tuple(sorted(pos['white']))
        initial_black_pos = tuple(sorted(pos['black']))
        
        print(f"--- Configuration {name} ---")

        # Condition 1: Equal number of knights
        if len(initial_white_pos) != len(initial_black_pos):
            print(f"Result: INVALID - Unequal number of knights ({len(initial_white_pos)} white, {len(initial_black_pos)} black).")
            print("-" * (len(name) + 20) + "\n")
            continue
            
        # Goal state is the swap of the initial state
        goal_white_pos = initial_black_pos
        goal_black_pos = initial_white_pos
        
        # State in queue: (white_positions, black_positions, turn)
        # Turn: 0 for White, 1 for Black. White starts.
        start_state = (initial_white_pos, initial_black_pos, 0)
        
        q = collections.deque([start_state])
        visited = {start_state}
        
        solution_found = False
        
        # BFS algorithm to find a path to the goal state
        while q:
            current_w_pos, current_b_pos, turn = q.popleft()
            
            # Check if current state is the goal state
            if current_w_pos == goal_white_pos and current_b_pos == goal_black_pos:
                solution_found = True
                break

            # Determine whose turn it is and get their pieces
            is_white_turn = turn == 0
            current_player_pos = current_w_pos if is_white_turn else current_b_pos
            other_player_pos = current_b_pos if is_white_turn else current_w_pos
            
            occupied_squares = set(current_w_pos) | set(current_b_pos)
            
            # Generate next possible states
            for i, p_start in enumerate(current_player_pos):
                for p_end in MOVES[p_start]:
                    if p_end not in occupied_squares:
                        # Create the new state by moving the piece
                        new_player_pos_list = list(current_player_pos)
                        new_player_pos_list[i] = p_end
                        new_player_pos = tuple(sorted(new_player_pos_list))
                        
                        next_w_pos = new_player_pos if is_white_turn else other_player_pos
                        next_b_pos = other_player_pos if is_white_turn else new_player_pos
                        
                        new_state = (next_w_pos, next_b_pos, 1 - turn)
                        
                        if new_state not in visited:
                            visited.add(new_state)
                            q.append(new_state)
                            
        if solution_found:
            print("Result: SOLVABLE")
            solvable_configs.append(name)
        else:
            print("Result: NOT Solvable")
        
        print("-" * (len(name) + 20) + "\n")

    print("Summary:")
    if solvable_configs:
        print(f"The solvable configurations are: {', '.join(sorted(solvable_configs))}")
    else:
        print("None of the configurations are solvable.")
    
solve_knights_puzzle()