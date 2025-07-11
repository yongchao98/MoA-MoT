import collections

def solve_knights_puzzle():
    """
    Solves the Knights Puzzle for the five given configurations on a 4x3 board.
    
    This script determines which initial setups can reach the goal state where
    white and black knights have swapped their starting positions.
    """

    # Step 1: Define the board and pre-compute knight moves
    # The 4x3 board is represented by squares numbered 0 to 11.
    # 0  1  2
    # 3  4  5
    # 6  7  8
    # 9 10 11
    
    def get_knight_moves():
        """Pre-computes all possible knight moves on the 4x3 board."""
        moves = collections.defaultdict(list)
        rows, cols = 4, 3
        for r in range(rows):
            for c in range(cols):
                pos = r * cols + c
                # All 8 possible L-shaped move offsets
                for dr, dc in [(1, 2), (1, -2), (-1, 2), (-1, -2),
                               (2, 1), (2, -1), (-2, 1), (-2, -1)]:
                    nr, nc = r + dr, c + dc
                    if 0 <= nr < rows and 0 <= nc < cols:
                        next_pos = nr * cols + nc
                        moves[pos].append(next_pos)
        return dict(moves)

    MOVES = get_knight_moves()

    # Step 2: Define the initial and goal configurations for A-E
    # A state is a tuple: (frozenset_of_white_knights, frozenset_of_black_knights)
    CONFIGS = {
        'A': {
            'initial': (frozenset({2, 5, 8, 11}), frozenset({0, 3, 6, 9})),
            'goal':    (frozenset({0, 3, 6, 9}), frozenset({2, 5, 8, 11})),
        },
        'B': {
            'initial': (frozenset({4, 9, 11}), frozenset({1, 5, 6})),
            'goal':    (frozenset({1, 5, 6}), frozenset({4, 9, 11})),
        },
        'C': {
            'initial': (frozenset({0, 7}), frozenset({2, 5})),
            'goal':    (frozenset({2, 5}), frozenset({0, 7})),
        },
        'D': {
            'initial': (frozenset({0, 7}), frozenset({4, 10})),
            'goal':    (frozenset({4, 10}), frozenset({0, 7})),
        },
        'E': {
            'initial': (frozenset({1, 2, 5}), frozenset({0, 3, 4})),
            'goal':    (frozenset({0, 3, 4}), frozenset({1, 2, 5})),
        }
    }

    def run_bfs_solver(initial_state, goal_state):
        """
        Solves the puzzle for a single configuration using Breadth-First Search.
        A search node is: (white_knights_pos, black_knights_pos, turn)
        where turn=0 for White and turn=1 for Black.
        """
        # White starts, so initial turn is 0
        start_node = (initial_state[0], initial_state[1], 0)
        
        queue = collections.deque([start_node])
        visited = {start_node}
        
        while queue:
            current_white, current_black, turn = queue.popleft()

            # Check if we've reached the goal board state
            if (current_white, current_black) == goal_state:
                return True

            # Determine whose turn it is and generate next moves
            knights_to_move = current_white if turn == 0 else current_black
            next_turn = 1 - turn
            occupied_squares = current_white.union(current_black)

            for start_pos in knights_to_move:
                for end_pos in MOVES.get(start_pos, []):
                    # A move is valid only to an empty square
                    if end_pos not in occupied_squares:
                        if turn == 0:  # White moved
                            new_white = (current_white - {start_pos}) | {end_pos}
                            new_black = current_black
                        else:  # Black moved
                            new_white = current_white
                            new_black = (current_black - {start_pos}) | {end_pos}
                        
                        next_node = (new_white, new_black, next_turn)

                        if next_node not in visited:
                            visited.add(next_node)
                            queue.append(next_node)
                            
        return False # Queue is empty, goal not found

    # Step 3: Iterate through each configuration and solve it
    solvable_configs = []
    print("Analyzing Knights Puzzle configurations...\n")
    for name, data in sorted(CONFIGS.items()):
        print(f"Checking configuration {name}...")
        is_solvable = run_bfs_solver(data['initial'], data['goal'])
        if is_solvable:
            solvable_configs.append(name)
        print(f"Result: Configuration {name} is {'solvable' if is_solvable else 'unsolvable'}.\n")
    
    # Step 4: Print the final summary
    print("="*35)
    print("           SUMMARY")
    print("="*35)
    if solvable_configs:
        print("The following configurations are solvable:")
        print(", ".join(solvable_configs))
    else:
        print("None of the configurations were found to be solvable.")
    print("="*35)

# Run the solver
solve_knights_puzzle()