import collections

def solve_sliding_puzzle():
    """
    Solves the sliding block puzzle to find the minimum moves to get the red piece to the top-left corner.

    The puzzle state is defined by the top-left coordinates of each of the 4 pieces.
    - R: Red 2x2 square
    - L: Large beige 2x3 rectangle
    - V: Vertical beige 2x1 rectangle
    - H: Horizontal beige 1x2 rectangle

    The board is a 4x4 grid. We use a Breadth-First Search (BFS) to find the shortest path.
    A move consists of sliding one piece one unit into an adjacent empty space.
    """
    
    # Piece dimensions (height, width)
    piece_dims = {
        'R': (2, 2),  # Red piece (target)
        'L': (2, 3),  # Large blocking piece
        'V': (2, 1),  # Vertical piece
        'H': (1, 2)   # Horizontal piece
    }

    # Initial state: positions of the top-left corner of each piece, derived from the puzzle image.
    # The state is represented as a tuple of piece positions in the order (R, L, V, H).
    initial_state = ((2, 2), (0, 0), (2, 0), (0, 3))
    
    # Goal state: Red piece R's top-left corner is at (0, 0)
    goal_r_pos = (0, 0)

    # BFS setup
    # The queue stores tuples of (state, number_of_moves)
    queue = collections.deque([(initial_state, 0)])
    # The visited set stores states we have already processed to avoid cycles and redundant work.
    visited = {initial_state}
    
    # Board dimensions
    board_height = 4
    board_width = 4

    while queue:
        current_state, moves = queue.popleft()
        
        # Unpack current positions from the state tuple
        r_pos, l_pos, v_pos, h_pos = current_state
        
        # Check if the current state is the goal state
        if r_pos == goal_r_pos:
            print(moves)
            return

        # --- Generate all valid next states from the current state ---
        
        # Map piece names to their current state details for easier iteration
        pieces = {
            'R': {'pos': r_pos, 'dim': piece_dims['R']},
            'L': {'pos': l_pos, 'dim': piece_dims['L']},
            'V': {'pos': v_pos, 'dim': piece_dims['V']},
            'H': {'pos': h_pos, 'dim': piece_dims['H']},
        }
        
        # Determine all occupied cells on the board to check for valid moves
        all_occupied = set()
        for piece_name, details in pieces.items():
            r, c = details['pos']
            h, w = details['dim']
            for i in range(h):
                for j in range(w):
                    all_occupied.add((r + i, c + j))

        # Try to move each of the 4 pieces
        for piece_idx, (name_to_move, details_to_move) in enumerate(pieces.items()):
            pos_to_move = details_to_move['pos']
            dim_to_move = details_to_move['dim']

            # Determine which cells are occupied by all *other* pieces
            current_piece_cells = set()
            r_start, c_start = pos_to_move
            h_piece, w_piece = dim_to_move
            for i in range(h_piece):
                for j in range(w_piece):
                    current_piece_cells.add((r_start + i, c_start + j))
            occupied_by_others = all_occupied - current_piece_cells

            # Try moving the selected piece one step in each of the 4 directions
            for dr, dc in [(0, 1), (0, -1), (1, 0), (-1, 0)]: # Right, Left, Down, Up
                new_r, new_c = pos_to_move[0] + dr, pos_to_move[1] + dc
                
                # Check if the new position is within the board boundaries
                if not (0 <= new_r and (new_r + dim_to_move[0]) <= board_height and \
                        0 <= new_c and (new_c + dim_to_move[1]) <= board_width):
                    continue
                
                # Determine the cells the piece would occupy in its new position
                new_piece_cells = set()
                for i in range(dim_to_move[0]):
                    for j in range(dim_to_move[1]):
                        new_piece_cells.add((new_r + i, new_c + j))

                # A move is valid if the new cells do not collide with any other piece
                if new_piece_cells.isdisjoint(occupied_by_others):
                    # If the move is valid, create the new state
                    new_positions = list(current_state)
                    new_positions[piece_idx] = (new_r, new_c)
                    new_state = tuple(new_positions)
                    
                    # If this new state has not been visited, add it to the queue and visited set
                    if new_state not in visited:
                        visited.add(new_state)
                        queue.append((new_state, moves + 1))

# Execute the solver function to find and print the minimum number of moves
solve_sliding_puzzle()