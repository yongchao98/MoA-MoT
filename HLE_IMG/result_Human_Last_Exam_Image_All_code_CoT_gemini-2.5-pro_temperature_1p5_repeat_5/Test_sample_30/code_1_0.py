import collections

# Board dimensions
ROWS, COLS = 4, 3

# Precompute all valid knight moves for each square to speed up the search.
# A dictionary where keys are (row, col) tuples and values are lists of valid move destinations.
VALID_MOVES = {}

def precompute_moves():
    """Calculates and stores all possible knight moves from each square on the board."""
    moves = [(-2, -1), (-2, 1), (-1, -2), (-1, 2),
             (1, -2), (1, 2), (2, -1), (2, 1)]
    for r in range(ROWS):
        for c in range(COLS):
            VALID_MOVES[(r, c)] = []
            for dr, dc in moves:
                nr, nc = r + dr, c + dc
                if 0 <= nr < ROWS and 0 <= nc < COLS:
                    VALID_MOVES[(r, c)].append((nr, nc))

class KnightsPuzzleSolver:
    """
    Solves the Knights Puzzle for a given initial configuration using Breadth-First Search.
    """
    def __init__(self, initial_config):
        """
        Initializes the solver with the starting positions of white and black knights.
        
        Args:
            initial_config (dict): A dictionary with keys 'W' and 'B', containing
                                   sets of (row, col) tuples for white and black knights.
        """
        # Use frozensets for hashability, required for the 'visited' set.
        self.initial_white = frozenset(initial_config['W'])
        self.initial_black = frozenset(initial_config['B'])
        
        # The goal is to swap the initial positions.
        self.goal_white = self.initial_black
        self.goal_black = self.initial_white

    def solve(self):
        """
        Performs a Breadth-First Search to find if the goal state is reachable.
        
        Returns:
            bool: True if the puzzle is solvable, False otherwise.
        """
        # A state is (white_positions, black_positions, turn_to_move)
        # White always moves first.
        initial_state = (self.initial_white, self.initial_black, 'W')
        
        # Queue for BFS, initialized with the starting state.
        q = collections.deque([initial_state])
        
        # Set to keep track of visited states to avoid cycles and redundant work.
        visited = {initial_state}
        
        while q:
            white_pos, black_pos, turn = q.popleft()
            
            # Check if the current configuration matches the goal.
            # The turn does not matter for the goal state itself, only the positions.
            if white_pos == self.goal_white and black_pos == self.goal_black:
                return True
                
            # Determine whose turn it is and who moves next.
            movers_pos = white_pos if turn == 'W' else black_pos
            next_turn = 'B' if turn == 'W' else 'W'
            
            occupied_squares = white_pos.union(black_pos)
            
            # Generate all possible next states from the current state.
            for start_pos in movers_pos:
                for end_pos in VALID_MOVES[start_pos]:
                    if end_pos not in occupied_squares:
                        # Create the new set of positions for the moved piece
                        new_movers_pos = (movers_pos - {start_pos}) | {end_pos}
                        
                        # Construct the next state tuple
                        if turn == 'W':
                            next_state = (frozenset(new_movers_pos), black_pos, next_turn)
                        else: # turn == 'B'
                            next_state = (white_pos, frozenset(new_movers_pos), next_turn)
                        
                        # If we haven't seen this state before, add it to the queue and visited set.
                        if next_state not in visited:
                            visited.add(next_state)
                            q.append(next_state)
                            
        # If the queue is empty and goal was not found, the puzzle is unsolvable.
        return False

def main():
    """
    Main function to define puzzle configurations, run the solver for each,
    and print the results.
    """
    # Precompute moves before starting the solvers.
    precompute_moves()

    # Define the five initial configurations from the image.
    # (row, col) coordinates are 0-indexed.
    CONFIGS = {
        'A': {
            'W': {(0, 2), (1, 2), (2, 2), (3, 2)},
            'B': {(0, 0), (1, 0), (2, 0), (3, 0)}
        },
        'B': {
            'W': {(1, 1), (3, 0), (3, 2)},
            'B': {(0, 0), (2, 0), (2, 1)}
        },
        'C': {
            'W': {(0, 1), (2, 1)},
            'B': {(0, 2), (1, 2)}
        },
        'D': {
            'W': {(0, 1), (2, 1)},
            'B': {(1, 1), (3, 0)}
        },
        'E': {
            'W': {(0, 1), (0, 2), (1, 2)},
            'B': {(0, 0), (1, 0), (1, 1)}
        }
    }
    
    solvable_configs = []
    print("Analyzing configurations...")
    for name, config in CONFIGS.items():
        solver = KnightsPuzzleSolver(config)
        is_solvable = solver.solve()
        print(f"Configuration {name}: {'Solvable' if is_solvable else 'Not Solvable'}")
        if is_solvable:
            solvable_configs.append(name)
    
    print("\nResult:")
    print(f"The solvable configurations are: {', '.join(sorted(solvable_configs))}")


if __name__ == '__main__':
    main()