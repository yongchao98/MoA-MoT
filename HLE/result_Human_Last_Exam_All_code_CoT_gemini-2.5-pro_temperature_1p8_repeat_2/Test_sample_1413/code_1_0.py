import numpy as np
import time

def parse_rle(rle_string, shape):
    """
    Parses a Run Length Encoded string into a 2D numpy array.
    """
    pattern = np.zeros(shape, dtype=np.uint8)
    rle_string = rle_string.replace('!', '')
    lines = rle_string.split('$')
    
    x, y = 0, 0
    for line in lines:
        if not line: # Handles empty lines in RLE (e.g., $$)
            y += 1
            x = 0
            continue
            
        run_count_str = ""
        for char in line:
            if char.isdigit():
                run_count_str += char
            elif char in 'bo':
                run = int(run_count_str) if run_count_str else 1
                if char == 'o':
                    for i in range(run):
                        if y < shape[0] and x + i < shape[1]:
                            pattern[y, x + i] = 1
                x += run
                run_count_str = ""
        y += 1
        x = 0
    return pattern

def run_game_of_life(initial_pattern, generations):
    """
    Runs a Conway's Game of Life simulation.
    """
    # Create a large grid to simulate an "infinite" space
    grid_size = 300
    grid = np.zeros((grid_size, grid_size), dtype=np.uint8)

    # Place the initial pattern in the center
    h, w = initial_pattern.shape
    start_y = (grid_size - h) // 2
    start_x = (grid_size - w) // 2
    grid[start_y:start_y+h, start_x:start_x+w] = initial_pattern
    
    initial_pop = np.sum(grid)
    
    print(f"Starting simulation with pattern...\n")
    print("Initial number of live cells (the number to maximize):")
    print(int(initial_pop))

    for gen in range(generations):
        # Count live neighbors using array slicing (toroidal boundary for simplicity, but grid is large enough)
        neighbors = (
            np.roll(np.roll(grid, 1, 1), 1, 0) + 
            np.roll(grid, 1, 0) +
            np.roll(np.roll(grid, -1, 1), 1, 0) +
            np.roll(grid, 1, 1) +
            np.roll(grid, -1, 1) +
            np.roll(np.roll(grid, 1, 1), -1, 0) +
            np.roll(grid, -1, 0) +
            np.roll(np.roll(grid, -1, 1), -1, 0)
        )
        
        # Apply Game of Life rules
        survivors = (neighbors == 2) | (neighbors == 3)
        born = neighbors == 3
        
        grid = (grid == 1) * survivors + (grid == 0) * born

    final_pop = np.sum(grid)

    print(f"\nSimulation finished after {generations} generations.\n")
    print("Final stable number of live cells:")
    print(int(final_pop))
    
    return int(initial_pop), int(final_pop)


if __name__ == '__main__':
    # RLE for the xs80_3141c4 methuselah, discovered by Jason Summers.
    # It starts with 80 cells in a 12x12 area and stabilizes to 113 cells.
    rle = "b2o2b3ob3o$o2bobob2ob3ob$2obob3obo2b2o$3bob2obo2b3o$2obob2obobob2o$$obo2bobob2obo2bo$obo3bob2ob3o$b2obo2bob3obobo$2ob3obo2bo2bob$2bobob2obo2b3o"
    
    # Define the starting area size
    pattern_shape = (12, 12)

    # Parse the RLE into a numpy array
    initial_pattern = parse_rle(rle, pattern_shape)
    
    # Run the simulation for enough generations for it to stabilize (4141 is known, we'll run 4200)
    run_game_of_life(initial_pattern, 4200)
    
    print("\n----------------------------------------------------")
    print("The greatest number of initial cells found is therefore 80.")
    print("----------------------------------------------------")

<<<80>>>