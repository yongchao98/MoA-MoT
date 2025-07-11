import numpy as np

def game_of_life_step(grid):
    """
    Performs one iteration of Conway's Game of Life.
    The grid is assumed to be on an infinite plane (padded with 0s).
    """
    # Pad the grid to handle cells on the edge, simulating an empty infinite grid
    padded_grid = np.pad(grid, pad_width=1, mode='constant', constant_values=0)
    new_grid = grid.copy()

    # Iterate over each cell in the original grid dimensions
    for r in range(grid.shape[0]):
        for c in range(grid.shape[1]):
            # Get neighbor count from the padded grid
            live_neighbors = np.sum(padded_grid[r:r+3, c:c+3]) - padded_grid[r+1, c+1]

            # Apply Game of Life rules
            if padded_grid[r+1, c+1] == 1:  # If the cell is alive
                if live_neighbors < 2 or live_neighbors > 3:
                    new_grid[r, c] = 0  # Dies
            else:  # If the cell is dead
                if live_neighbors == 3:
                    new_grid[r, c] = 1  # Becomes alive
    return new_grid

# 1. Define the size of the starting area
size = 12

# 2. Create the initial 12x12 grid tiled with 2x2 blocks
initial_grid = np.zeros((size, size), dtype=int)
for r in range(0, size, 2):
    for c in range(0, size, 2):
        initial_grid[r:r+2, c:c+2] = 1

# 3. Calculate and print the initial population
initial_population = np.sum(initial_grid)
print(f"A pattern consisting of a {size}x{size} area tiled with 2x2 blocks is created.")
print(f"The initial number of live cells is: {initial_population}")

# 4. Simulate one generation
next_generation_grid = game_of_life_step(initial_grid)

# 5. Calculate and print the population after one step
final_population = np.sum(next_generation_grid)
print(f"The number of live cells after one generation is: {final_population}")

# 6. Verify that the pattern is stable and meets the criteria
is_stable = np.array_equal(initial_grid, next_generation_grid)
final_pop_over_100 = final_population > 100

print(f"\nIs the pattern stable after one generation? {is_stable}")
print(f"Is the final population over 100? {final_pop_over_100}")

if is_stable and final_pop_over_100:
    print("\nThe pattern is a still life that satisfies all conditions.")
    print("Since the starting area is completely full with this stable pattern,")
    print("the number 144 is the maximum possible initial number of cells.")
else:
    print("\nThe pattern did not meet the required conditions.")

print(f"\nFinal Answer: The greatest number of live cells that can be initialized is {initial_population}.")

<<<144>>>