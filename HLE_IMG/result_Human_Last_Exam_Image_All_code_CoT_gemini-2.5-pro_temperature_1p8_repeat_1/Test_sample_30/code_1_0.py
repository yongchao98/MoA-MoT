import collections

def solve_knights_puzzle():
    """
    This script solves the Knights Puzzle for five different initial configurations
    on a 4x3 chessboard. It determines which of these configurations are solvable
    by performing a Breadth-First Search of the state space.
    """

    # --- Step 1: Define the game board and rules ---
    # The 4x3 board is represented by squares indexed 0 to 11.
    # (row, col) -> row * 3 + col
    # We pre-calculate all possible knight moves from each square.
    MOVES = {}
    ROWS, COLS = 4, 3
    for r in range(ROWS):
        for c in range(COLS):
            idx = r * COLS + c
            MOVES[idx] = []
            # All 8 possible L-shaped moves for a knight
            for dr, dc in [(1, 2), (1, -2), (-1, 2), (-1, -2),
                           (2, 1), (2, -1), (-2, 1), (-2, -1)]:
                nr, nc = r + dr, c + dc
                if 0 <= nr < ROWS and 0 <= nc < COLS:
                    n_idx = nr * COLS + nc
                    MOVES[idx].append(n_idx)
    MOVES = {k: tuple(v) for k, v in MOVES.items()}

    # --- Step 2: Define the BFS solver ---
    def is_solvable(initial_white, initial_black):
        """
        Performs a Breadth-First Search (BFS) to find if the goal state is reachable.
        A state is defined by the positions of all knights and whose turn it is.
        """
        # Use sorted tuples for a canonical representation of positions
        initial_white = tuple(sorted(initial_white))
        initial_black = tuple(sorted(initial_black))

        initial_state = (initial_white, initial_black)
        goal_state = (initial_black, initial_white)
        
        # Trivial case: solvable in 0 moves if start is goal.
        if initial_state == goal_state:
            return True

        # The queue stores tuples of (state, turn)
        queue = collections.deque([(initial_state, 'W')])
        # The visited set stores (state, turn) to avoid cycles
        visited = {(initial_state, 'W')}

        while queue:
            current_state, turn = queue.popleft()

            if current_state == goal_state:
                return True

            current_white, current_black = current_state
            occupied_squares = set(current_white) | set(current_black)

            pieces_to_move, next_turn = (current_white, 'B') if turn == 'W' else (current_black, 'W')
            
            # Generate all possible next states
            for i, piece_pos in enumerate(pieces_to_move):
                for dest_pos in MOVES[piece_pos]:
                    if dest_pos not in occupied_squares:
                        # Create the new configuration after the move
                        new_pieces_list = list(pieces_to_move)
                        new_pieces_list[i] = dest_pos
                        new_pieces_tuple = tuple(sorted(new_pieces_list))
                        
                        next_state = (new_pieces_tuple, current_black) if turn == 'W' else (current_white, new_pieces_tuple)

                        if (next_state, next_turn) not in visited:
                            visited.add((next_state, next_turn))
                            queue.append((next_state, next_turn))
        
        # If the queue is empty and goal was not found, it's unsolvable.
        return False

    # --- Step 3: Define the initial configurations from the image ---
    configs = {
        'A': {
            'W': (2, 5, 8, 11),  # White at (0,2), (1,2), (2,2), (3,2)
            'B': (0, 3, 6, 9),   # Black at (0,0), (1,0), (2,0), (3,0)
        },
        'B': {
            'W': (4, 9, 11),     # White at (1,1), (3,0), (3,2)
            'B': (1, 6, 8),      # Black at (0,1), (2,0), (2,2)
        },
        'C': {
            'W': (0, 7),         # White at (0,0), (2,1)
            'B': (2, 5),         # Black at (0,2), (1,2)
        },
        'D': {
            'W': (1, 7),         # White at (0,1), (2,1)
            'B': (4, 10),        # Black at (1,1), (3,1)
        },
        'E': {
            'W': (1, 2, 5),      # White at (0,1), (0,2), (1,2)
            'B': (0, 3, 4),      # Black at (0,0), (1,0), (1,1)
        }
    }
    
    # --- Step 4: Solve for each configuration and print results ---
    print("Analyzing the Knights Puzzle configurations...")
    solvable_configs = []
    for name, pos in configs.items():
        result = is_solvable(pos['W'], pos['B'])
        print(f"Configuration {name}: {'Solvable' if result else 'Unsolvable'}")
        if result:
            solvable_configs.append(name)
            
    final_answer = ", ".join(sorted(solvable_configs))
    print(f"\nThe solvable configurations are: {final_answer}")

# Execute the solver
solve_knights_puzzle()