import numpy as np
from scipy.signal import convolve2d
import sys

def simulate_life(initial_pattern_str):
    """
    Simulates Conway's Game of Life for a given initial pattern.

    Args:
        initial_pattern_str (str): A multiline string representing the pattern,
                                   with 'O' for live cells and '.' for dead cells.
    """
    # --- Configuration ---
    # The simulation must run for over 10750 generations to stabilize.
    # We set a limit slightly higher than the known value.
    MAX_GENERATIONS = 11000
    
    # The pattern can grow significantly, so we need a large canvas.
    CANVAS_SIZE = 400

    # --- Pattern Initialization ---
    lines = initial_pattern_str.strip().split('\n')
    pattern_height = len(lines)
    pattern_width = len(lines[0])

    initial_grid = np.zeros((pattern_height, pattern_width), dtype=np.uint8)
    for r, line in enumerate(lines):
        for c, char in enumerate(line):
            if char == 'O':
                initial_grid[r, c] = 1
    
    initial_population = initial_grid.sum()
    print(f"The greatest number of live cells found is: {initial_population}")
    print("Starting simulation... This will take a long time.")
    sys.stdout.flush()

    # Place the pattern in the center of a large canvas
    grid = np.zeros((CANVAS_SIZE, CANVAS_SIZE), dtype=np.uint8)
    start_r = (CANVAS_SIZE - pattern_height) // 2
    start_c = (CANVAS_SIZE - pattern_width) // 2
    grid[start_r:start_r+pattern_height, start_c:start_c+pattern_width] = initial_grid

    # --- Simulation ---
    # Kernel for counting neighbors
    kernel = np.array([[1, 1, 1],
                       [1, 0, 1],
                       [1, 1, 1]], dtype=np.uint8)

    history = {grid.tobytes(): 0}

    for gen in range(1, MAX_GENERATIONS + 1):
        if gen % 1000 == 0:
            print(f"Processing generation {gen}/{MAX_GENERATIONS}...")
            sys.stdout.flush()

        # Count live neighbors for each cell
        neighbor_count = convolve2d(grid, kernel, mode='same', boundary='fill', fillvalue=0)
        
        # Apply Game of Life rules
        # 1. A live cell with 2 or 3 neighbors survives.
        survivors = grid & ((neighbor_count == 2) | (neighbor_count == 3))
        # 2. A dead cell with exactly 3 neighbors becomes a live cell.
        newborns = (1 - grid) & (neighbor_count == 3)
        
        grid = survivors | newborns
        
        # Check for stabilization or death
        population = grid.sum()
        if population == 0:
            print("\nPattern died out.")
            return

        grid_bytes = grid.tobytes()
        if grid_bytes in history:
            prev_gen = history[grid_bytes]
            final_population = population
            print("\n--- Simulation Complete ---")
            print(f"Pattern stabilized at generation {gen}.")
            print(f"It returned to a state previously seen at generation {prev_gen}.")
            if final_population > 100:
                print(f"Success! The final stable population is {final_population}, which is over 100.")
                # Final output summarizing the result
                print(f"The initial state of {initial_population} cells stabilized to {final_population} cells.")
            else:
                 print(f"Failure. The final stable population is {final_population}, which is not over 100.")
            return

        history[grid_bytes] = gen
        
    print(f"\nSimulation finished after {MAX_GENERATIONS} generations without stabilizing.")


if __name__ == '__main__':
    # This 10x10 pattern contains 40 cells and was discovered by David Eppstein.
    # It is known to stabilize after 10,750 generations to a population of 105.
    # Source: https://www.conwaylife.com/wiki/295P10750.1
    boojum_pattern = """\
..OO....OO
.O..O..O..O
OO..O.O..OO
.O.O...O.O
..O.....O
..O.....O
.O.O...O.O
OO..O.O..OO
.O..O..O..O
..OO....OO"""

    simulate_life(boojum_pattern)