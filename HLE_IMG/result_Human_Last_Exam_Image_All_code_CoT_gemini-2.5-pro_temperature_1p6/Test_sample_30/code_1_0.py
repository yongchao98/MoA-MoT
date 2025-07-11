import collections

def solve_knights_puzzle():
    """
    Solves the Knights Puzzle for the five given configurations.
    
    The solution uses a multi-step process:
    1.  A "color parity" check to find configurations that are impossible due
        to the nature of knight moves on a colored board.
    2.  An "initial stalemate" check to find configurations where the first
        player has no moves.
    3.  A Breadth-First Search (BFS) for the remaining configurations to
        definitively prove solvability by finding a solution path.
    """
    # ---- Configuration Data ----
    # Board: 4x3 grid. Coords (row, col) from (0,0) to (3,2).
    CONFIGS = {
        'A': {'B': frozenset([(0, 0), (1, 0), (2, 0), (3, 0)]), 'W': frozenset([(0, 2), (1, 2), (2, 2), (3, 2)])},
        'B': {'B': frozenset([(0, 1), (2, 0), (2, 2)]), 'W': frozenset([(1, 1), (3, 0), (3, 2)])},
        'C': {'B': frozenset([(0, 2), (1, 2)]), 'W': frozenset([(0, 0), (2, 1)])},
        'D': {'B': frozenset([(1, 1), (3, 1)]), 'W': frozenset([(0, 0), (2, 0)])},
        'E': {'B': frozenset([(0, 0), (1, 0), (1, 1)]), 'W': frozenset([(0, 1), (0, 2), (1, 2)])}
    }

    # ---- Helper Functions ----
    def get_knight_moves_map(rows, cols):
        moves_map = collections.defaultdict(list)
        knight_hops = [(1, 2), (1, -2), (-1, 2), (-1, -2), (2, 1), (2, -1), (-2, 1), (-2, -1)]
        for r in range(rows):
            for c in range(cols):
                for dr, dc in knight_hops:
                    nr, nc = r + dr, c + dc
                    if 0 <= nr < rows and 0 <= nc < cols:
                        moves_map[(r, c)].append((nr, nc))
        return moves_map

    KNIGHT_MOVES = get_knight_moves_map(4, 3)

    def get_board_state(black_pos, white_pos):
        return (tuple(sorted(list(black_pos))), tuple(sorted(list(white_pos))))

    # ---- Puzzle Solvability Checks ----
    def check_color_parity(initial_pos, target_pos):
        def count_on_dark_squares(positions):
            return sum(1 for r, c in positions if (r + c) % 2 == 0)

        initial_b_dark = count_on_dark_squares(initial_pos['B'])
        target_b_dark = count_on_dark_squares(target_pos['B'])
        b_moves_odd = (initial_b_dark % 2) != (target_b_dark % 2)

        initial_w_dark = count_on_dark_squares(initial_pos['W'])
        target_w_dark = count_on_dark_squares(target_pos['W'])
        w_moves_odd = (initial_w_dark % 2) != (target_w_dark % 2)

        if w_moves_odd and b_moves_odd:
            return False, "Both colors require an odd number of moves, which is impossible."
        
        return True, "Passes color parity check."

    def check_stalemate(initial_pos):
        white_knights = initial_pos['W']
        occupied_squares = initial_pos['B'].union(white_knights)
        for knight_pos in white_knights:
            for move_pos in KNIGHT_MOVES[knight_pos]:
                if move_pos not in occupied_squares:
                    return False, "White has valid moves."
        return True, "White has no valid moves (stalemate)."

    def bfs_solve(initial_pos, target_pos):
        initial_state = (get_board_state(initial_pos['B'], initial_pos['W']), 'W')
        target_state_board = get_board_state(target_pos['B'], target_pos['W'])
        queue = collections.deque([initial_state])
        visited = {initial_state}
        
        # Limit search depth to prevent excessive runtime, though not strictly needed for these small problems
        max_depth = 20
        depth_map = {initial_state: 0}

        while queue:
            current_board_tuple, turn = queue.popleft()
            depth = depth_map[(current_board_tuple, turn)]
            
            if depth >= max_depth: continue

            if current_board_tuple == target_state_board:
                return True, f"Solution found via BFS at depth {depth}."

            current_black_pos = frozenset(current_board_tuple[0])
            current_white_pos = frozenset(current_board_tuple[1])
            occupied = current_black_pos.union(current_white_pos)
            
            player_pos = current_white_pos if turn == 'W' else current_black_pos
            next_turn = 'B' if turn == 'W' else 'W'
            
            for knight_pos in player_pos:
                for move_pos in KNIGHT_MOVES[knight_pos]:
                    if move_pos not in occupied:
                        new_player_pos = player_pos - {knight_pos} | {move_pos}
                        next_b = current_black_pos if turn == 'W' else new_player_pos
                        next_w = new_player_pos if turn == 'W' else current_white_pos
                        next_state = (get_board_state(next_b, next_w), next_turn)
                        if next_state not in visited:
                            visited.add(next_state)
                            depth_map[next_state] = depth + 1
                            queue.append(next_state)
                            
        return False, "No solution exists (proven by exhaustive BFS)."

    # Main analysis loop
    solvable_configs = []
    for name in sorted(CONFIGS.keys()):
        is_solvable = True
        reason = ""
        initial_pos = CONFIGS[name]
        target_pos = {'W': initial_pos['B'], 'B': initial_pos['W']}
        
        parity_ok, r = check_color_parity(initial_pos, target_pos)
        if not parity_ok:
            is_solvable, reason = False, r
        else:
            stalemated, r = check_stalemate(initial_pos)
            if stalemated:
                is_solvable, reason = False, r
            else:
                bfs_ok, r = bfs_solve(initial_pos, target_pos)
                if not bfs_ok:
                    is_solvable, reason = False, r
                else:
                    reason = r

        if is_solvable:
            solvable_configs.append(name)
    
    print("The following configurations are solvable:")
    print(", ".join(solvable_configs))


solve_knights_puzzle()
<<<A, D>>>