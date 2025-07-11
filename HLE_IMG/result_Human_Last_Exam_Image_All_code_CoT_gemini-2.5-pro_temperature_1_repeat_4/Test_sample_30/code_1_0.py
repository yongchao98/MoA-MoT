import collections

# Global constants and caches for board dimensions and knight moves
BOARD_ROWS = 4
BOARD_COLS = 3
_knight_moves_cache = {}

def get_knight_moves(r, c):
    """Generates all valid knight moves from (r, c) on a 4x3 board."""
    if (r, c) in _knight_moves_cache:
        return _knight_moves_cache[(r, c)]

    moves = []
    deltas = [(-2, -1), (-2, 1), (-1, -2), (-1, 2),
              (1, -2), (1, 2), (2, -1), (2, 1)]
    for dr, dc in deltas:
        nr, nc = r + dr, c + dc
        if 0 <= nr < BOARD_ROWS and 0 <= nc < BOARD_COLS:
            moves.append((nr, nc))
    _knight_moves_cache[(r, c)] = moves
    return moves

def solve_puzzle(name, config):
    """
    Solves the knight puzzle for a given configuration using BFS.
    Prints whether the configuration is solvable or not.
    """
    initial_white = config['W']
    initial_black = config['B']

    # Condition 1: Must be same number of knights to swap all positions.
    if len(initial_white) != len(initial_black):
        print(f"Configuration {name}: Unsolvable (has {len(initial_white)} White vs {len(initial_black)} Black knights).")
        return False

    # The goal is to swap the initial positions.
    goal_white = frozenset(initial_black)
    goal_black = frozenset(initial_white)
    
    start_white = frozenset(initial_white)
    start_black = frozenset(initial_black)
    
    # State in queue: (white_positions, black_positions, turn)
    # Positions are frozensets to be hashable for the visited set.
    queue = collections.deque([(start_white, start_black, 'W')])  # White moves first
    visited = {(start_white, start_black, 'W')}

    # Condition 2: Can the first player (White) make any move?
    can_white_move = False
    occupied_at_start = start_white.union(start_black)
    for r_w, c_w in start_white:
        for nr, nc in get_knight_moves(r_w, c_w):
            if (nr, nc) not in occupied_at_start:
                can_white_move = True
                break
        if can_white_move:
            break
    if not can_white_move:
        print(f"Configuration {name}: Unsolvable (White has no opening move).")
        return False

    # Set a reasonable search limit to prevent very long runtimes.
    max_states_to_check = 2_000_000 
    count = 0

    while queue:
        count += 1
        if count > max_states_to_check:
            print(f"Configuration {name}: Unsolvable (search space is too large or no solution exists).")
            return False

        current_w_pos, current_b_pos, turn = queue.popleft()

        # Check for goal state
        if current_w_pos == goal_white and current_b_pos == goal_black:
            print(f"Configuration {name}: Solvable (found a path of moves).")
            return True

        occupied = current_w_pos.union(current_b_pos)
        
        movers, next_turn = (current_w_pos, 'B') if turn == 'W' else (current_b_pos, 'W')

        # Iterate through each knight of the current player
        for r, c in movers:
            # Iterate through all its possible moves
            for nr, nc in get_knight_moves(r, c):
                # A move is valid only if the destination square is empty
                if (nr, nc) not in occupied:
                    if turn == 'W':
                        new_w_pos = (current_w_pos - {(r, c)}) | {(nr, nc)}
                        new_b_pos = current_b_pos
                    else: # turn == 'B'
                        new_w_pos = current_w_pos
                        new_b_pos = (current_b_pos - {(r, c)}) | {(nr, nc)}
                    
                    new_state_tuple = (new_w_pos, new_b_pos, next_turn)

                    if new_state_tuple not in visited:
                        visited.add(new_state_tuple)
                        queue.append((new_w_pos, new_b_pos, next_turn))
                        
    print(f"Configuration {name}: Unsolvable (no path to goal state).")
    return False

def main():
    """
    Defines the puzzle configurations from the image and solves them.
    """
    # Coordinates are (row, col) with (0,0) at top-left.
    configs = {
        'A': {
            'B': {(0, 0), (1, 0), (2, 0), (3, 0)},
            'W': {(0, 2), (1, 2), (2, 2), (3, 2)}
        },
        'B': {
            'B': {(0, 1), (1, 0), (2, 0), (2, 2)},
            'W': {(1, 1), (3, 0), (3, 2)}
        },
        'C': {
            'B': {(0, 2), (1, 2)},
            'W': {(0, 0), (2, 1)}
        },
        'D': {
            'B': {(1, 1), (3, 1)},
            'W': {(0, 0), (2, 1)}
        },
        'E': {
            'B': {(0, 0), (1, 0), (1, 1)},
            'W': {(0, 1), (0, 2), (1, 2)}
        }
    }
    
    print("Analyzing Knights Puzzle configurations...\n")
    solvable_configs = []
    for name, config_data in sorted(configs.items()):
        if solve_puzzle(name, config_data):
            solvable_configs.append(name)
    
    print("\n--- Summary ---")
    if solvable_configs:
        # The final answer format requires printing the letters.
        final_answer = ", ".join(sorted(solvable_configs))
        print(f"The solvable configurations are: {final_answer}")
    else:
        print("None of the configurations were found to be solvable.")

if __name__ == '__main__':
    main()