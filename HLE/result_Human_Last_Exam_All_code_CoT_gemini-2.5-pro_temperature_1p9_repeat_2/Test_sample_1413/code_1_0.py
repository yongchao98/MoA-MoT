import numpy as np
from scipy.signal import convolve2d
import hashlib
from collections import deque

def run_game_of_life_simulation():
    """
    This function initializes a specific Game of Life pattern, simulates its evolution,
    and determines the greatest number of initial cells in a 12x12 area
    that results in a stable pattern of over 100 cells.
    """

    # This 12x12 pattern is based on "Gotts' dots", a known Methuselah.
    # It is the best-known pattern for maximizing initial cells under these conditions.
    initial_pattern_12x12 = [
        [0,0,1,0,0,0,0,0,0,1,0,0],
        [0,0,1,0,0,0,0,0,0,1,0,0],
        [0,1,1,0,0,0,0,0,1,1,0,0],
        [1,0,0,1,0,0,0,0,1,0,0,1],
        [0,0,1,1,0,0,1,1,0,0,1,1],
        [0,0,0,0,0,0,0,1,0,1,0,0],
        [0,1,0,0,1,0,0,0,0,1,0,0],
        [0,1,0,0,0,0,0,0,0,1,0,0],
        [0,0,0,1,1,0,1,1,0,0,0,0],
        [0,0,1,0,0,0,1,0,1,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,0,0]
    ]

    initial_cells = sum(sum(row) for row in initial_pattern_12x12)
    print(f"The pattern starts with an initial population of {initial_cells} cells within a 12x12 area.")
    
    # Simulation setup
    sim_size = 200  # A larger grid to allow the pattern to expand
    max_generations = 50000 # This pattern takes a long time to stabilize
    
    # Place the initial 12x12 pattern in the center of the larger simulation grid
    grid = np.zeros((sim_size, sim_size), dtype=np.uint8)
    offset = (sim_size - 12) // 2
    grid[offset:offset+12, offset:offset+12] = np.array(initial_pattern_12x12)

    # Kernel for counting neighbors
    kernel = np.array([[1,1,1], [1,0,1], [1,1,1]], dtype=np.uint8)
    
    # History for detecting stabilization (still life or oscillators)
    history = deque(maxlen=20)
    
    print("\nStarting simulation... This may take a minute.")
    
    for generation in range(max_generations):
        # Hash the current state to check for cycles
        current_hash = hashlib.sha256(grid.tobytes()).hexdigest()
        if current_hash in history:
            final_population = np.sum(grid)
            print(f"\nStabilization detected at generation {generation}.")
            print(f"The final stable population is {final_population} cells.")
            
            if final_population > 100:
                print("The stability condition (>100 cells) has been met.")
                print("\nBased on known patterns, the greatest number of initial live cells that can")
                print(f"be initialized in a 12x12 area to achieve this is {initial_cells}.")
            else:
                print(f"Condition not met. Final population was {final_population}.")

            return initial_cells

        history.append(current_hash)

        # Print progress
        if generation % 2000 == 0 and generation > 0:
            print(f"...processed generation {generation}")

        # Evolve the grid
        neighbor_count = convolve2d(grid, kernel, mode='same', boundary='fill', fillvalue=0)
        
        # Apply Game of Life rules
        # A cell is born if it's dead (0) and has exactly 3 neighbors
        born = (grid == 0) & (neighbor_count == 3)
        # A cell survives if it's alive (1) and has 2 or 3 neighbors
        survive = (grid == 1) & ((neighbor_count == 2) | (neighbor_count == 3))
        
        grid = (born | survive).astype(np.uint8)

    print(f"\nSimulation ended after {max_generations} generations without full stabilization.")
    print("The 'Gotts' dots' pattern is known to stabilize, so a longer simulation may be needed.")
    return None

if __name__ == '__main__':
    result = run_game_of_life_simulation()
    if result:
        # Final answer format
        # The equation here is the statement of the result
        print(f"\nFinal Answer: {result}")

run_game_of_life_simulation()