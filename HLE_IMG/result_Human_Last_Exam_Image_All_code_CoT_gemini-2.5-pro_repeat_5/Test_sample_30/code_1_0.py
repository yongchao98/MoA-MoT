import collections

def solve_knights_puzzle():
    """
    Solves the Knights Puzzle for five different configurations.

    The function defines the initial and goal states for five puzzle configurations (A-E)
    on 4x3 or 3x4 chessboards. It uses a Breadth-First Search (BFS) algorithm to
    determine if the goal state (positions of white and black knights swapped) is
    reachable from the initial state, following the rules of alternating moves
    (white starts).

    The state in the search is represented by a tuple containing the board layout
    and the player whose turn it is. The board itself is a tuple where each element
    represents a square: 1 for a white knight, -1 for a black knight, and 0 for empty.

    The function iterates through each configuration, runs the BFS solver, and
    collects the labels of the solvable puzzles. Finally, it prints the result.
    """

    def generate_adj(rows, cols):
        """Generates an adjacency list for knight moves on a board of given size."""
        adj = collections.defaultdict(list)
        moves = [(1, 2), (1, -2), (-1, 2), (-1, -2),
                 (2, 1), (2, -1), (-2, 1), (-2, -1)]
        for r in range(rows):
            for c in range(cols):
                pos = r * cols + c
                for dr, dc in moves:
                    nr, nc = r + dr, c + dc
                    if 0 <= nr < rows and 0 <= nc < cols:
                        npos = nr * cols + nc
                        adj[pos].append(npos)
        return adj

    def bfs_solver(config):
        """
        Performs a Breadth-First Search to find if a configuration is solvable.
        """
        initial_board = config['initial']
        goal_board = config['goal']
        rows = config['rows']
        cols = config['cols']
        
        adj = generate_adj(rows, cols)
        
        # State: (board_tuple, turn), where turn is 1 for white, -1 for black.
        # White always starts.
        start_state = (initial_board, 1)
        
        queue = collections.deque([start_state])
        visited = {start_state}
        
        # Set a reasonable limit for the number of states to explore.
        max_states = 200000
        count = 0

        while queue:
            count += 1
            if count > max_states:
                return False  # Assumed unsolvable if search is too large

            current_board, current_turn = queue.popleft()
            
            if current_board == goal_board:
                return True
                
            # Find all knights of the color whose turn it is
            knight_indices = [i for i, piece in enumerate(current_board) if piece == current_turn]
            
            for pos in knight_indices:
                for new_pos in adj[pos]:
                    # A knight can only move to an empty square
                    if current_board[new_pos] == 0:
                        next_board_list = list(current_board)
                        next_board_list[new_pos] = current_turn
                        next_board_list[pos] = 0
                        next_board = tuple(next_board_list)
                        
                        next_turn = -current_turn
                        next_state = (next_board, next_turn)
                        
                        if next_state not in visited:
                            visited.add(next_state)
                            queue.append(next_state)
                            
        return False

    # --- Define the 5 configurations ---
    # Representation: 1=White, -1=Black, 0=Empty
    # Board squares are numbered 0-11, row by row.

    # Config A: 4x3 board
    initial_A = tuple([-1, 0, 1] * 4)
    goal_A = tuple([1, 0, -1] * 4)

    # Config B: 4x3 board
    # Black at (0,0), (2,0) -> pos 0, 6
    # White at (1,1), (3,1) -> pos 4, 10
    b_list = [0]*12; b_list[0]=-1; b_list[6]=-1; b_list[4]=1; b_list[10]=1
    initial_B = tuple(b_list)
    b_goal_list = [0]*12; b_goal_list[0]=1; b_goal_list[6]=1; b_goal_list[4]=-1; b_goal_list[10]=-1
    goal_B = tuple(b_goal_list)

    # Config C: 4x3 board
    # Black at (0,2), (1,2) -> pos 2, 5
    # White at (0,0), (2,1) -> pos 0, 7
    c_list = [0]*12; c_list[2]=-1; c_list[5]=-1; c_list[0]=1; c_list[7]=1
    initial_C = tuple(c_list)
    c_goal_list = [0]*12; c_goal_list[2]=1; c_goal_list[5]=1; c_goal_list[0]=-1; c_goal_list[7]=-1
    goal_C = tuple(c_goal_list)
    
    # Config D: 4x3 board
    # Black at (1,1), (3,1) -> pos 4, 10
    # White at (0,0), (2,0) -> pos 0, 6
    d_list = [0]*12; d_list[4]=-1; d_list[10]=-1; d_list[0]=1; d_list[6]=1
    initial_D = tuple(d_list)
    d_goal_list = [0]*12; d_goal_list[4]=1; d_goal_list[10]=1; d_goal_list[0]=-1; d_goal_list[6]=-1
    goal_D = tuple(d_goal_list)

    # Config E: 3x4 board
    # Black at (0,0),(0,1),(1,0),(1,1) -> pos 0,1,4,5
    # White at (0,2),(0,3),(1,2),(1,3) -> pos 2,3,6,7
    e_list = [0]*12; e_list[0]=-1;e_list[1]=-1;e_list[4]=-1;e_list[5]=-1
    e_list[2]=1;e_list[3]=1;e_list[6]=1;e_list[7]=1
    initial_E = tuple(e_list)
    e_goal_list = [0]*12; e_goal_list[0]=1;e_goal_list[1]=1;e_goal_list[4]=1;e_goal_list[5]=1
    e_goal_list[2]=-1;e_goal_list[3]=-1;e_goal_list[6]=-1;e_goal_list[7]=-1
    goal_E = tuple(e_goal_list)

    configurations = {
        'A': {'initial': initial_A, 'goal': goal_A, 'rows': 4, 'cols': 3},
        'B': {'initial': initial_B, 'goal': goal_B, 'rows': 4, 'cols': 3},
        'C': {'initial': initial_C, 'goal': goal_C, 'rows': 4, 'cols': 3},
        'D': {'initial': initial_D, 'goal': goal_D, 'rows': 4, 'cols': 3},
        'E': {'initial': initial_E, 'goal': goal_E, 'rows': 3, 'cols': 4}
    }

    solvable_configs = []
    for name, config_data in configurations.items():
        if bfs_solver(config_data):
            solvable_configs.append(name)
    
    print("The solvable configurations are:")
    print(', '.join(sorted(solvable_configs)))

solve_knights_puzzle()