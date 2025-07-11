import collections

def solve_knights_puzzle():
    """
    Solves the Knights Puzzle for 5 initial configurations on a 4x3 board.
    The function defines the configurations, runs a Breadth-First Search for each,
    and prints which ones are solvable.
    """

    # Helper function to get valid knight moves from a given position
    def get_knight_moves(pos, rows=4, cols=3):
        r, c = pos
        moves = []
        deltas = [(1, 2), (1, -2), (-1, 2), (-1, -2),
                  (2, 1), (2, -1), (-2, 1), (-2, -1)]
        for dr, dc in deltas:
            nr, nc = r + dr, c + dc
            if 0 <= nr < rows and 0 <= nc < cols:
                moves.append((nr, nc))
        return moves

    # The main BFS solver for a single configuration
    def is_solvable(initial_white, initial_black):
        goal_white = frozenset(initial_black)
        goal_black = frozenset(initial_white)

        start_state = (frozenset(initial_white), frozenset(initial_black), 'W')
        goal_positions = (goal_white, goal_black)

        queue = collections.deque([start_state])
        visited = {start_state}
        
        # Using a generous iteration limit as a safeguard
        max_iterations = 250000 
        count = 0

        while queue and count < max_iterations:
            count += 1
            current_white_pos, current_black_pos, turn = queue.popleft()

            if (current_white_pos, current_black_pos) == goal_positions:
                return True

            occupied_squares = current_white_pos.union(current_black_pos)
            next_turn = 'B' if turn == 'W' else 'W'
            
            knights_to_move = current_white_pos if turn == 'W' else current_black_pos

            for knight_pos in knights_to_move:
                for move_pos in get_knight_moves(knight_pos):
                    if move_pos not in occupied_squares:
                        new_white_set = set(current_white_pos)
                        new_black_set = set(current_black_pos)
                        
                        if turn == 'W':
                            new_white_set.remove(knight_pos)
                            new_white_set.add(move_pos)
                        else:
                            new_black_set.remove(knight_pos)
                            new_black_set.add(move_pos)
                        
                        new_state = (frozenset(new_white_set), frozenset(new_black_set), next_turn)

                        if new_state not in visited:
                            visited.add(new_state)
                            queue.append(new_state)
        return False

    # Define the 5 configurations based on the image (row, col) from top-left (0,0)
    configs = {
        'A': {'black': {(0, 0), (1, 0), (2, 0), (3, 0)}, 'white': {(0, 2), (1, 2), (2, 2), (3, 2)}},
        'B': {'black': {(0, 1), (1, 0), (2, 0), (2, 2)}, 'white': {(1, 1), (3, 0), (3, 1), (3, 2)}},
        'C': {'black': {(0, 2), (1, 2)}, 'white': {(0, 0), (2, 1)}},
        'D': {'black': {(1, 1), (3, 1)}, 'white': {(0, 0), (2, 1)}},
        'E': {'black': {(0, 0), (1, 0), (1, 1)}, 'white': {(0, 1), (0, 2), (1, 2)}}
    }

    solvable_configs = []
    print("Analyzing the Knight's Puzzle configurations...")
    for name in sorted(configs.keys()):
        config = configs[name]
        if is_solvable(config['white'], config['black']):
            solvable_configs.append(name)
            print(f"Configuration {name} is solvable.")
        else:
            print(f"Configuration {name} is unsolvable.")
    
    print("\nSummary:")
    if solvable_configs:
        print("The following configurations are solvable: " + ", ".join(solvable_configs))
    else:
        print("None of the configurations are solvable.")

# Execute the solver
solve_knights_puzzle()