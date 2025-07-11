def solve_game_of_life_puzzle():
    """
    Calculates how many of the 512 possible 3x3 starting configurations in
    Conway's Game of Life eventually die out on an infinite plane.
    """

    def get_next_generation(live_cells):
        """
        Calculates the next state of the grid based on the rules of the game.

        Args:
            live_cells: A set of (row, col) tuples representing live cells.

        Returns:
            A new set of (row, col) tuples for the next generation.
        """
        # We only need to check cells that are currently alive or are neighbors of living cells.
        potential_cells = set()
        for r, c in live_cells:
            for dr in range(-1, 2):
                for dc in range(-1, 2):
                    potential_cells.add((r + dr, c + dc))

        new_live_cells = set()
        for r, c in potential_cells:
            # Count live neighbors for each potential cell
            live_neighbors = 0
            for dr in range(-1, 2):
                for dc in range(-1, 2):
                    if dr == 0 and dc == 0:
                        continue
                    if (r + dr, c + dc) in live_cells:
                        live_neighbors += 1
            
            # Apply the rules of Life: B3/S23
            is_currently_alive = (r, c) in live_cells
            if not is_currently_alive and live_neighbors == 3:
                # Rule for reproduction
                new_live_cells.add((r, c))
            elif is_currently_alive and (live_neighbors == 2 or live_neighbors == 3):
                # Rule for survival
                new_live_cells.add((r, c))
                
        return new_live_cells

    total_configs = 2**9
    dies_out_count = 0
    # A safe generation limit to prevent extremely long-running (but not infinite) simulations
    # from slowing down the process. Patterns from a 3x3 grid stabilize or die out quickly.
    max_generations = 150 

    # Iterate through all 2^9 = 512 possible configurations
    for i in range(total_configs):
        # Generate the initial configuration from the integer 'i'
        # The j-th bit of 'i' corresponds to the state of the j-th cell in the 3x3 grid.
        live_cells = set()
        for j in range(9):
            if (i >> j) & 1:
                row = j // 3
                col = j % 3
                live_cells.add((row, col))

        # Use history to detect cycles (oscillators, still lifes)
        history = {frozenset(live_cells)}
        
        configuration_dies = False
        for _ in range(max_generations):
            if not live_cells:
                # The configuration has no more living cells
                configuration_dies = True
                break

            live_cells = get_next_generation(live_cells)
            
            # Use a frozenset (an immutable set) for the key in our history tracking
            frozen_state = frozenset(live_cells)
            if frozen_state in history:
                # A cycle is detected, so this pattern will never die out
                break
            history.add(frozen_state)
        
        if configuration_dies:
            dies_out_count += 1
            
    # Output the final equation as requested
    print(f"Total possible initial configurations for a 3x3 grid: {total_configs}")
    non_dying_configs = total_configs - dies_out_count
    print(f"Number of configurations that result in no living cells: {dies_out_count}")
    print(f"Number of configurations that result in stable or oscillating patterns: {non_dying_configs}")
    print(f"\nFinal Equation: {total_configs} (total) - {non_dying_configs} (stable/oscillating) = {dies_out_count} (die out)")


solve_game_of_life_puzzle()
<<<376>>>