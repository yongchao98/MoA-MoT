import numpy as np
from scipy.signal import convolve2d
import time

def run_game_of_life_and_find_stable_state(initial_pattern, max_generations=2000, stability_check=20):
    """
    Simulates Conway's Game of Life for a given initial pattern.

    Args:
        initial_pattern (list of tuples): A list of (row, col) coordinates for the initial live cells.
        max_generations (int): The maximum number of generations to simulate.
        stability_check (int): The number of consecutive generations with the same population
                                to consider the pattern stable.

    Returns:
        A tuple containing (initial_cell_count, final_cell_count, generations_to_stabilize).
        Returns (initial_cell_count, -1, -1) if it does not stabilize within max_generations.
    """
    initial_cell_count = len(initial_pattern)
    
    # Determine the required grid size from the pattern's bounding box.
    # We add a large buffer to allow for expansion.
    if not initial_pattern:
        return 0, 0, 0
    max_r = max(p[0] for p in initial_pattern)
    max_c = max(p[1] for p in initial_pattern)
    grid_size = max(max_r, max_c) + 300  # Large buffer for evolution
    
    grid = np.zeros((grid_size, grid_size), dtype=np.int8)

    # Place the pattern in the center of the grid
    offset_r = (grid_size - max_r) // 2
    offset_c = (grid_size - max_c) // 2
    for r, c in initial_pattern:
        grid[r + offset_r, c + offset_c] = 1

    # Kernel for counting neighbors
    kernel = np.array([[1, 1, 1],
                       [1, 0, 1],
                       [1, 1, 1]], dtype=np.int8)

    population_history = []

    print(f"Starting simulation with {initial_cell_count} cells...")
    
    for generation in range(max_generations):
        current_population = np.sum(grid)
        population_history.append(current_population)

        # Check for stability
        if generation > stability_check:
            last_n_counts = population_history[-stability_check:]
            if len(set(last_n_counts)) == 1:
                print(f"Pattern stabilized at generation {generation - stability_check + 1}.")
                final_cell_count = int(current_population)
                return initial_cell_count, final_cell_count, generation - stability_check + 1
        
        # Calculate neighbors using 2D convolution
        neighbors = convolve2d(grid, kernel, mode='same', boundary='fill', fillvalue=0)
        
        # Apply the rules of Life
        survivors = (grid == 1) & ((neighbors == 2) | (neighbors == 3))
        births = (grid == 0) & (neighbors == 3)
        
        grid = (survivors | births).astype(np.int8)

    print("Simulation finished without reaching a stable state within max_generations.")
    return initial_cell_count, -1, -1

if __name__ == "__main__":
    # The "h-to-r converter" pattern.
    # It has 19 cells and fits within an 11x6 bounding box.
    # It evolves into an R-pentomino, which stabilizes at 116 cells.
    h_to_r_converter_pattern = [
        (0, 8), (0, 9),
        (1, 7), (1, 9),
        (2, 2), (2, 3), (2, 7), (2, 9),
        (3, 1), (3, 3), (3, 6), (3, 10),
        (4, 1), (4, 3), (4, 7), (4, 9),
        (5, 0), (5, 3), (5, 8)
    ]

    start_time = time.time()
    initial_cells, final_cells, generations = run_game_of_life_and_find_stable_state(
        h_to_r_converter_pattern,
        max_generations=1500  # The pattern stabilizes around generation 1200
    )
    end_time = time.time()

    print("\n--- Results ---")
    if final_cells > 100:
        print(f"The chosen pattern meets the criteria.")
        print(f"Initial number of live cells: {initial_cells}")
        print(f"Final stable number of live cells: {final_cells}")
        print(f"This is the greatest known number of initial cells for a pattern in a 12x12 area that stabilizes to over 100 cells.")
    elif final_cells != -1:
        print(f"The chosen pattern stabilized at {final_cells} cells, which does not meet the >100 condition.")
    else:
        print("The simulation did not find a stable state.")
    
    print(f"\n(Simulation took {end_time - start_time:.2f} seconds)")
