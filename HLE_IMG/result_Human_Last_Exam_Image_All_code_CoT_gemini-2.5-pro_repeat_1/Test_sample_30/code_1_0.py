import collections

def solve_knights_puzzle():
    """
    Solves the Knights Puzzle for the 5 given configurations on a 4x3 board.
    It uses a theoretical parity argument to quickly rule out impossible cases,
    and then uses Breadth-First Search (BFS) to confirm the solvability of the remaining cases.
    """

    # Board layout:
    # 0  1  2
    # 3  4  5
    # 6  7  8
    # 9 10 11
    ADJ = {
        0: [5, 7], 1: [6, 8], 2: [3, 7],
        3: [2, 8, 10], 4: [9, 11], 5: [0, 10],
        6: [1, 11], 7: [0, 2], 8: [1, 3],
        9: [4, 8], 10: [3, 5], 11: [4, 6]
    }

    # Configurations from the image
    # W: white knight positions, B: black knight positions
    CONFIGS = {
        'A': {'W': {2, 5, 8, 11}, 'B': {0, 3, 6, 9}},
        'B': {'W': {3, 5, 10}, 'B': {1, 6, 8}},
        'C': {'W': {0, 4}, 'B': {1, 5}},
        'D': {'W': {0, 7, 8}, 'B': {1, 4, 9}},
        'E': {'W': {1, 2, 5}, 'B': {0, 3, 4}}
    }
    
    solvable_configs = []

    print("Analyzing configurations...")
    for label, config in CONFIGS.items():
        num_knights_per_color = len(config['W'])
        
        # Parity argument: The number of knight pairs (N) must be even.
        if num_knights_per_color % 2 != 0:
            print(f"Configuration {label}: Unsolvable (has {num_knights_per_color} pairs of knights, which is an odd number).")
            continue

        # Use BFS to find a solution for cases with an even number of knight pairs.
        initial_white = frozenset(config['W'])
        initial_black = frozenset(config['B'])
        
        goal_white = initial_black
        goal_black = initial_white

        # State: (current_white_positions, current_black_positions, turn)
        # turn=1 for white, turn=-1 for black
        start_state = (initial_white, initial_black, 1) # White starts

        queue = collections.deque([(start_state, 0)])
        visited = {start_state}
        
        # Limit search depth to keep the search practical.
        max_depth = 16 
        found = False

        while queue:
            (current_white, current_black, turn), depth = queue.popleft()
            
            if depth >= max_depth:
                continue

            # Check for goal state. A solution path must have an even number of moves.
            if current_white == goal_white and current_black == goal_black and depth % 2 == 0:
                print(f"Configuration {label}: Solvable in {depth} moves.")
                solvable_configs.append(label)
                found = True
                break

            # Generate next valid states
            next_turn = -turn
            movable_knights = current_white if turn == 1 else current_black
            static_knights = current_black if turn == 1 else current_white
            occupied_squares = movable_knights.union(static_knights)

            for pos in movable_knights:
                for move_to in ADJ[pos]:
                    if move_to not in occupied_squares:
                        new_movable = movable_knights.difference({pos}).union({move_to})
                        
                        if turn == 1:
                            new_state = (new_movable, static_knights, next_turn)
                        else:
                            new_state = (static_knights, new_movable, next_turn)

                        if new_state not in visited:
                            visited.add(new_state)
                            queue.append((new_state, depth + 1))
        
        if not found:
            print(f"Configuration {label}: Not solved within {max_depth} moves.")

    print("\nBased on the analysis, the solvable configurations are:")
    print(", ".join(sorted(solvable_configs)))

solve_knights_puzzle()