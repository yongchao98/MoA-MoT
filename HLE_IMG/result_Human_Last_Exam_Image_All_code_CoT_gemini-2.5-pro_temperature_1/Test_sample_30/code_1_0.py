import collections

def solve_knights_puzzle():
    """
    Solves the Knights Puzzle for five initial configurations on a 4x3 board.
    """

    # --- Setup ---
    ROWS, COLS = 4, 3
    MOVES = [(1, 2), (1, -2), (-1, 2), (-1, -2),
             (2, 1), (2, -1), (-2, 1), (-2, -1)]
    
    # Board indices:
    # 0  1  2
    # 3  4  5
    # 6  7  8
    # 9 10 11
    # Use 1 for White, -1 for Black, 0 for Empty

    configs = {
        "A": (-1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1),
        "B": (-1, 0, 0, 0, 1, 0, -1, 0, -1, 1, 0, 1),
        "C": (0, 1, -1, 0, 0, -1, 0, 1, 0, 0, 0, 0),
        "D": (0, 1, 0, 0, -1, 0, 0, 1, 0, 0, -1, 0),
        "E": (-1, 1, 1, -1, -1, 1, 0, 0, 0, 0, 0, 0)
    }

    # --- Helper Functions ---
    def index_to_coords(i):
        return (i // COLS, i % COLS)

    def coords_to_index(r, c):
        return r * COLS + c

    def get_goal_state(initial_state):
        return tuple(-p for p in initial_state)

    def get_next_states(board, player):
        next_boards = []
        board_list = list(board)
        for i, piece in enumerate(board_list):
            if piece == player:
                r, c = index_to_coords(i)
                for dr, dc in MOVES:
                    nr, nc = r + dr, c + dc
                    if 0 <= nr < ROWS and 0 <= nc < COLS:
                        ni = coords_to_index(nr, nc)
                        if board_list[ni] == 0:
                            new_board_list = board_list[:]
                            new_board_list[i], new_board_list[ni] = new_board_list[ni], new_board_list[i]
                            next_boards.append(tuple(new_board_list))
        return next_boards

    # --- Solver Function ---
    def is_solvable(name, initial_state):
        print(f"--- Analyzing Configuration {name} ---")
        
        # 1. Graph Connectivity Check
        c1_indices = {0, 1, 2, 3, 4, 5, 7, 9, 11}
        
        white_in_c1 = sum(1 for i, p in enumerate(initial_state) if p == 1 and i in c1_indices)
        black_in_c1 = sum(1 for i, p in enumerate(initial_state) if p == -1 and i in c1_indices)
        
        num_white = initial_state.count(1)
        num_black = initial_state.count(-1)
        
        white_in_c2 = num_white - white_in_c1
        black_in_c2 = num_black - black_in_c1

        if white_in_c1 != black_in_c1 or white_in_c2 != black_in_c2:
            print(f"Result: Unsolvable. The number of white and black knights in disconnected graph components do not match.")
            print(f"Component 1: {white_in_c1} White vs {black_in_c1} Black")
            print(f"Component 2: {white_in_c2} White vs {black_in_c2} Black")
            return False

        print("Connectivity check passed. Starting state-space search...")
        
        # 2. Breadth-First Search (BFS)
        goal_state = get_goal_state(initial_state)
        queue = collections.deque([(initial_state, 1)]) # (board, player_to_move)
        visited = {(initial_state, 1)}

        max_iterations = 200000 # Safety break
        count = 0
        
        while queue:
            if count > max_iterations:
                print(f"Result: Search timed out after {max_iterations} states. Assuming unsolvable.")
                return False
            count += 1
            
            current_board, player_to_move = queue.popleft()

            if current_board == goal_state:
                print(f"Result: Solvable. Goal state found after checking {count} states.")
                return True

            next_player = -player_to_move
            for next_board in get_next_states(current_board, player_to_move):
                if (next_board, next_player) not in visited:
                    visited.add((next_board, next_player))
                    queue.append((next_board, next_player))

        print(f"Result: Unsolvable. Search completed after checking {count} states, no solution found.")
        return False

    # --- Main Execution ---
    solvable_configs = []
    for name, initial_state in configs.items():
        if is_solvable(name, initial_state):
            solvable_configs.append(name)
        print("")

    print("=========================================")
    print(f"Final Answer: The solvable configurations are: {', '.join(solvable_configs)}")
    print("<<<" + ", ".join(solvable_configs) + ">>>")

solve_knights_puzzle()