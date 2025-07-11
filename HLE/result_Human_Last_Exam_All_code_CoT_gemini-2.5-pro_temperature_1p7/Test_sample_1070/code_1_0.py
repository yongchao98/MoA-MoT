import collections

def solve_chess_endgame():
    """
    Solves the described chess endgame problem using Breadth-First Search (BFS)
    to find the fastest checkmate, considering the special rule for the black king.
    """

    # --- Helper Functions ---

    def to_coords(alg_notation):
        """Converts algebraic notation (e.g., 'd5') to 0-indexed coordinates (e.g., (3, 4))."""
        file = ord(alg_notation[0]) - ord('a')
        rank = int(alg_notation[1]) - 1
        return (file, rank)

    def to_alg(coords):
        """Converts 0-indexed coordinates (e.g., (3, 4)) to algebraic notation (e.g., 'd5')."""
        file = chr(coords[0] + ord('a'))
        rank = str(coords[1] + 1)
        return file + rank

    def is_on_board(pos):
        """Checks if a position (file, rank) is on the 8x8 board."""
        return 0 <= pos[0] <= 7 and 0 <= pos[1] <= 7

    def is_white_square(pos):
        """
        Checks if a square is white.
        a1=(0,0) is black (sum=0). Squares where (file_index + rank_index) is odd are white.
        """
        return (pos[0] + pos[1]) % 2 == 1

    # --- Move Generation and State Analysis ---

    def get_attacked_squares(wk_pos, wb_pos, bk_pos):
        """
        Returns a set of all squares attacked by white's king and bishop.
        The bishop's line of sight is blocked by the other pieces.
        """
        attacked = set()
        # White King attacks
        for dx in [-1, 0, 1]:
            for dy in [-1, 0, 1]:
                if dx == 0 and dy == 0:
                    continue
                pos = (wk_pos[0] + dx, wk_pos[1] + dy)
                if is_on_board(pos):
                    attacked.add(pos)

        # White Bishop attacks
        for dx in [-1, 1]:
            for dy in [-1, 1]:
                for i in range(1, 8):
                    pos = (wb_pos[0] + i * dx, wb_pos[1] + i * dy)
                    if not is_on_board(pos):
                        break
                    attacked.add(pos)
                    # The bishop's attack is blocked by any piece on its path.
                    if pos == wk_pos or pos == bk_pos:
                        break
        return attacked

    def get_white_moves(state):
        """Generates all possible next states from white's moves."""
        wk_pos, wb_pos, bk_pos = state
        moves = []
        
        # White King moves
        for dx in [-1, 0, 1]:
            for dy in [-1, 0, 1]:
                if dx == 0 and dy == 0:
                    continue
                new_wk_pos = (wk_pos[0] + dx, wk_pos[1] + dy)
                if not is_on_board(new_wk_pos):
                    continue
                if new_wk_pos == wb_pos:  # Cannot move to square of own piece
                    continue
                # King cannot move adjacent to the opponent's king
                if abs(new_wk_pos[0] - bk_pos[0]) <= 1 and abs(new_wk_pos[1] - bk_pos[1]) <= 1:
                    continue
                move_str = f"K{to_alg(new_wk_pos)}"
                moves.append(((new_wk_pos, wb_pos, bk_pos), move_str))

        # White Bishop moves
        for dx in [-1, 1]:
            for dy in [-1, 1]:
                for i in range(1, 8):
                    new_wb_pos = (wb_pos[0] + i * dx, wb_pos[1] + i * dy)
                    if not is_on_board(new_wb_pos):
                        break
                    # The bishop is blocked by any piece on its path.
                    if new_wb_pos == wk_pos or new_wb_pos == bk_pos:
                        break
                    move_str = f"B{to_alg(new_wb_pos)}"
                    moves.append(((wk_pos, new_wb_pos, bk_pos), move_str))
        return moves

    def get_black_moves(state):
        """Generates all possible next states from black's (cursed) moves."""
        wk_pos, wb_pos, bk_pos = state
        attacked_by_white = get_attacked_squares(wk_pos, wb_pos, bk_pos)
        moves = []

        # Black King moves (to adjacent, un-attacked, WHITE squares)
        for dx in [-1, 0, 1]:
            for dy in [-1, 0, 1]:
                if dx == 0 and dy == 0:
                    continue
                new_bk_pos = (bk_pos[0] + dx, bk_pos[1] + dy)
                if not is_on_board(new_bk_pos):
                    continue
                if not is_white_square(new_bk_pos):  # The curse applies here
                    continue
                if new_bk_pos in attacked_by_white:
                    continue
                # move_str is not needed here but kept for consistency
                moves.append(((wk_pos, wb_pos, new_bk_pos), f"K{to_alg(new_bk_pos)}"))
        return moves

    # --- BFS Implementation ---
    
    initial_state = (to_coords('d2'), to_coords('e2'), to_coords('d5'))
    # The queue stores tuples of (state, path_of_moves)
    q = collections.deque([(initial_state, [])])
    # The visited set stores tuples of (state, turn) to handle repeated positions
    # Turn: 0 for White to move, 1 for Black to move
    visited = {(initial_state, 0)}

    while q:
        current_state, path = q.popleft()
        wk_pos, wb_pos, bk_pos = current_state
        turn = len(path) % 2

        if turn == 0:  # White's turn
            for new_state, move_str in get_white_moves(current_state):
                if (new_state, 1) not in visited:
                    new_path = path + [move_str]
                    q.append((new_state, new_path))
                    visited.add((new_state, 1))
        else:  # Black's turn
            is_in_check = bk_pos in get_attacked_squares(wk_pos, wb_pos, bk_pos)
            black_moves = get_black_moves(current_state)

            if is_in_check and not black_moves:
                # Checkmate found!
                first_move = path[0]
                # Number of moves for white is half the total ply, rounded up
                num_moves = (len(path) + 1) // 2
                
                # Print each part of the final answer as requested
                # The "equation" instruction is interpreted as printing the components.
                print(f"{first_move}, {num_moves}")
                return

            # If not a checkmate, continue the search
            for new_state, move_str in black_moves:
                if (new_state, 0) not in visited:
                    new_path = path + [move_str]
                    q.append((new_state, new_path))
                    visited.add((new_state, 0))

# Execute the solver
solve_chess_endgame()