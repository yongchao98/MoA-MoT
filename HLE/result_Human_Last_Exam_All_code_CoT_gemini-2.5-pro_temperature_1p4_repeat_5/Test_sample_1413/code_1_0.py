import numpy as np
import random
from collections import deque
import hashlib
import time

# --- Configuration ---
# The size of the starting area for the pattern.
START_AREA_SIZE = 12
# The size of the simulation grid. Must be large enough to contain pattern expansion.
# Known complex patterns can grow quite large.
SIMULATION_GRID_SIZE = 400
# The maximum number of generations to simulate before giving up on a pattern.
# Famous methuselahs can run for 10s of thousands of generations.
# Set high for a thorough search, lower for a quick test (e.g., 5000).
MAX_GENERATIONS = 100000
# Number of random patterns to try for each initial cell count.
# Higher values increase the chance of finding a rare pattern.
# Set to 10-50 for a quick test.
TRIALS_PER_COUNT = 100
# Population target for the stable pattern.
STABLE_POP_TARGET = 100
# How many past generations to store to check for stabilization (detects oscillators).
STABILIZATION_HISTORY_SIZE = 100


def create_random_pattern(width, height, num_cells):
    """Creates a random pattern with a specific number of live cells."""
    if num_cells > width * height:
        raise ValueError("Number of cells exceeds the area size.")
    
    pattern = np.zeros((height, width), dtype=np.uint8)
    # Get a list of all possible cell coordinates
    all_coords = [(r, c) for r in range(height) for c in range(width)]
    # Choose a random sample of coordinates to set to 1
    live_cell_coords = random.sample(all_coords, num_cells)
    for r, c in live_cell_coords:
        pattern[r, c] = 1
    return pattern

def step(grid):
    """Computes one generation of Conway's Game of Life."""
    # Pad the grid to handle boundary conditions by wrapping around (torus)
    padded_grid = np.pad(grid, 1, mode='constant')
    
    # Count live neighbors for each cell
    neighbor_count = np.zeros_like(grid)
    for r in range(grid.shape[0]):
        for c in range(grid.shape[1]):
            # Sum the 3x3 window in the padded grid, subtracting the cell itself
            neighbor_count[r, c] = padded_grid[r:r+3, c:c+3].sum() - grid[r, c]

    # Apply Game of Life rules
    # 1. A live cell with 2 or 3 neighbors survives.
    rule1 = (grid == 1) & ((neighbor_count == 2) | (neighbor_count == 3))
    # 2. A dead cell with exactly 3 neighbors becomes a live cell.
    rule2 = (grid == 0) & (neighbor_count == 3)
    
    return (rule1 | rule2).astype(np.uint8)

def hash_grid(grid):
    """Generates a hash for the current grid state."""
    return hashlib.sha256(grid.tobytes()).hexdigest()

def run_simulation(initial_pattern):
    """
    Runs a simulation for a given pattern.
    Returns (is_stable, final_pop, generations, final_pattern)
    """
    # Embed the initial pattern in the center of the larger simulation grid
    grid = np.zeros((SIMULATION_GRID_SIZE, SIMULATION_GRID_SIZE), dtype=np.uint8)
    start_r = (SIMULATION_GRID_SIZE - initial_pattern.shape[0]) // 2
    start_c = (SIMULATION_GRID_SIZE - initial_pattern.shape[1]) // 2
    grid[start_r:start_r+initial_pattern.shape[0], start_c:start_c+initial_pattern.shape[1]] = initial_pattern

    history = deque(maxlen=STABILIZATION_HISTORY_SIZE)
    
    for generation in range(MAX_GENERATIONS):
        # Optimization: if population is zero or too large, stop
        pop_count = np.sum(grid)
        if pop_count == 0:
            return True, 0, generation, grid

        # Check for stabilization
        h = hash_grid(grid)
        if h in history:
            final_pop = np.sum(grid)
            return True, final_pop, generation, grid
        history.append(h)
        
        grid = step(grid)

    # If loop finishes, it did not stabilize in time
    return False, np.sum(grid), MAX_GENERATIONS, grid

def search_for_methuselah():
    """Main search function."""
    max_initial_cells = START_AREA_SIZE * START_AREA_SIZE
    
    print("Starting search for the greatest number of initial live cells in a 12x12 area")
    print(f"that stabilizes to a population over {STABLE_POP_TARGET}.")
    print(f"This may take a very long time...\n")

    for initial_cells in range(max_initial_cells, 0, -1):
        print(f"--- Searching for a working pattern with {initial_cells} initial cells ---")
        for trial in range(TRIALS_PER_COUNT):
            if (trial + 1) % 10 == 0:
                print(f"  Trial {trial + 1}/{TRIALS_PER_COUNT}...")
                
            pattern = create_random_pattern(START_AREA_SIZE, START_AREA_SIZE, initial_cells)
            
            is_stable, final_pop, generations, final_grid = run_simulation(pattern)
            
            if is_stable and final_pop > STABLE_POP_TARGET:
                print("\n" + "="*50)
                print("!!! SUCCESS: Found a qualifying pattern !!!")
                print(f"Greatest initial cell count found: {initial_cells}")
                print(f"Number of generations to stabilize: {generations}")
                print(f"Final stable population: {final_pop}")
                print("Initial 12x12 Pattern ('1' is a live cell):")
                # Print the pattern row by row
                for row in pattern:
                    print(" ".join(map(str, row)))
                print("="*50 + "\n")
                return initial_cells
                
    print("\nSearch complete. No qualifying pattern was found with the given parameters.")
    return 0

if __name__ == '__main__':
    start_time = time.time()
    greatest_num = search_for_methuselah()
    end_time = time.time()
    
    if greatest_num > 0:
        print(f"The greatest number of initial cells found by this search is: {greatest_num}")
    
    print(f"Total search time: {end_time - start_time:.2f} seconds.")
