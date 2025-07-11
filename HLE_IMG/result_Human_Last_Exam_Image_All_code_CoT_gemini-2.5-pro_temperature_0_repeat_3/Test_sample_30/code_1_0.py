import collections

def solve_knights_puzzle():
    """
    Solves the Knights Puzzle for five given configurations on a 4x3 board
    by performing a Breadth-First Search (BFS) on the state space.
    """

    # --- 1. Define Board and Moves ---
    # The 4x3 board is indexed from 0 to 11, row by row.
    # 0  1  2
    # 3  4  5
    # 6  7  8
    # 9 10 11
    rows, cols = 4, 3
    moves = collections.defaultdict(list)
    for r in range(rows):
        for c in range(cols):
            pos = r * cols + c
            # All 8 possible knight move offsets
            for dr, dc in [(1, 2), (1, -2), (-1, 2), (-1, -2),
                           (2, 1), (2, -1), (-2, 1), (-2, -1)]:
                nr, nc = r + dr, c + dc
                if 0 <= nr < rows and 0 <= nc < cols:
                    npos = nr * cols + nc
                    moves[pos].append(npos)

    # --- 2. Define Initial Configurations ---
    # 'B': Black Knight, 'W': White Knight, 'E': Empty Square
    configs = {
        'A': ('B', 'E', 'W', 'B', 'E', 'W', 'B', 'E', 'W', 'B', 'E', 'W'),
        'B': ('E', 'B', 'E', 'W', 'B', 'E', 'B', 'W', 'B', 'W', 'E', 'W'),
        'C': ('W', 'E', 'B', 'E', 'E', 'B', 'E', 'W', 'E', 'E', 'E', 'E'),
        'D': ('W', 'E', 'E', 'E', 'B', 'E', 'E', 'W', 'E', 'E', 'B', 'E'),
        'E': ('B', 'W', 'W', 'B', 'B', 'W', 'E', 'E', 'E', 'E', 'E', 'E')
    }

    # --- 3. The Solver Function (BFS) ---
    def is_solvable(initial_state):
        """
        Determines if a given configuration is solvable using BFS.
        """
        # Derive the goal state by swapping knight colors
        goal_state_list = []
        for piece in initial_state:
            if piece == 'W':
                goal_state_list.append('B')
            elif piece == 'B':
                goal_state_list.append('W')
            else:
                goal_state_list.append('E')
        goal_state = tuple(goal_state_list)

        # The queue stores tuples of (board_state, current_turn)
        # White moves first
        queue = collections.deque([(initial_state, 'W')])
        
        # The visited set also stores (board_state, current_turn)
        visited = {(initial_state, 'W')}

        while queue:
            current_board, current_turn = queue.popleft()

            if current_board == goal_state:
                return True

            next_turn = 'B' if current_turn == 'W' else 'W'

            # Find all knights of the current color and generate next moves
            for pos, piece in enumerate(current_board):
                if piece == current_turn:
                    for dest in moves[pos]:
                        if current_board[dest] == 'E':
                            next_board_list = list(current_board)
                            next_board_list[pos], next_board_list[dest] = 'E', current_turn
                            next_board = tuple(next_board_list)

                            if (next_board, next_turn) not in visited:
                                visited.add((next_board, next_turn))
                                queue.append((next_board, next_turn))
        
        return False

    # --- 4. Solve for each configuration and print results ---
    solvable_configs = []
    print("Analyzing the five knight puzzle configurations...")
    for name, initial_state in sorted(configs.items()):
        is_sol = is_solvable(initial_state)
        print(f"Configuration {name}: {'Solvable' if is_sol else 'Not Solvable'}")
        if is_sol:
            solvable_configs.append(name)
    
    print("\nConclusion:")
    if not solvable_configs:
        print("None of the configurations are solvable.")
    else:
        print(f"The solvable configuration(s) are: {', '.join(solvable_configs)}")

# Execute the solver
solve_knights_puzzle()