import heapq

def solve_water_well():
    """
    Solves the water well problem by simulating the filling process.
    """
    grid = [
        [1, 5, 27, 22, 28, 40, 14],
        [39, 13, 17, 30, 41, 12, 2],
        [32, 35, 24, 25, 19, 47, 34],
        [16, 33, 10, 42, 7, 44, 18],
        [3, 8, 45, 37, 4, 21, 20],
        [15, 46, 38, 6, 26, 48, 49],
        [9, 23, 31, 29, 11, 36, 43],
    ]

    rows, cols = len(grid), len(grid[0])
    
    # Create a mapping from depth value to coordinates
    coords = {grid[r][c]: (r, c) for r in range(rows) for c in range(cols)}

    start_coords = coords[1]
    target_coords = coords[43]

    # Priority queue stores (depth, r, c) for cells on the boundary of the pool
    pq = []
    
    # Set of (r, c) tuples for cells that are part of the pool
    pool = set()
    
    total_time = 0
    current_level = 0
    
    # Start the process with the first cell
    heapq.heappush(pq, (grid[start_coords[0]][start_coords[1]], start_coords[0], start_coords[1]))

    print("Starting the simulation...")
    print("-" * 50)

    while pq:
        wall_height, r, c = heapq.heappop(pq)

        if (r, c) in pool:
            continue

        prev_pool_size = len(pool)
        
        # Calculate time to fill
        time_to_add = 0
        
        # Special case for the very first cell
        if not pool:
            time_to_add = wall_height
            print(f"Time to fill the first cell (depth {wall_height}): {time_to_add} minutes.")
            current_level = wall_height
        # For all subsequent cells
        elif wall_height > current_level:
            time_to_add = (wall_height - current_level) * prev_pool_size
            print(f"Time to raise level for {prev_pool_size} cells from {current_level}' to {wall_height}': ({wall_height} - {current_level}) * {prev_pool_size} = {time_to_add} minutes.")
            current_level = wall_height
        else: # wall_height <= current_level
            print(f"Water flows into cell with depth {wall_height}' (current level is {current_level}'). No extra time needed.")

        total_time += time_to_add
        print(f"New cell added: depth {grid[r][c]} at ({r},{c}). Pool size: {prev_pool_size + 1}. Current level: {current_level}'. Total time: {total_time:.0f} minutes.")
        print("-" * 50)

        # Add the new cell to the pool
        pool.add((r, c))

        # Check if the target cell has been reached
        if (r, c) == target_coords:
            print("\nWater has started to fill the target section 43.")
            break
            
        # Add new neighbors to the boundary
        for dr, dc in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
            nr, nc = r + dr, c + dc
            if 0 <= nr < rows and 0 <= nc < cols and (nr, nc) not in pool:
                heapq.heappush(pq, (grid[nr][nc], nr, nc))

    print(f"\nThe final answer is the total time elapsed until water begins to rise on section 43.")
    print(f"Total time = {total_time:.0f}")
    
# To run the simulation, uncomment the line below
solve_water_well()
print("<<<661>>>")