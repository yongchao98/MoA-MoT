import heapq

def solve_water_well():
    """
    Simulates the filling of the well to determine when the water level
    at section 43 begins to rise.
    """
    depths = [
        [1, 5, 27, 22, 28, 40, 14],
        [39, 13, 17, 30, 41, 12, 2],
        [32, 35, 24, 25, 19, 47, 34],
        [16, 33, 10, 42, 7, 44, 18],
        [3, 8, 45, 37, 4, 21, 20],
        [15, 46, 38, 6, 26, 48, 49],
        [9, 23, 31, 29, 11, 36, 43]
    ]
    rows, cols = 7, 7
    start_node = (0, 0)
    target_node = (6, 6)

    visited = [[False for _ in range(cols)] for _ in range(rows)]
    pq = []  # Min-heap: (wall_height, r, c)

    # Initial state
    pool_cells = set()
    water_level = depths[start_node[0]][start_node[1]]
    total_volume = 0
    
    volume_increments = []

    def add_to_pool_and_pq(r, c):
        """Helper to add a cell to the pool and its neighbors to the PQ."""
        if visited[r][c]:
            return
        
        visited[r][c] = True
        pool_cells.add((r,c))
        
        for dr, dc in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
            nr, nc = r + dr, c + dc
            if 0 <= nr < rows and 0 <= nc < cols and not visited[nr][nc]:
                wall_height = max(depths[r][c], depths[nr][nc])
                heapq.heappush(pq, (wall_height, nr, nc))

    # Start with the source node
    add_to_pool_and_pq(start_node[0], start_node[1])

    while pq:
        spill_level, r, c = heapq.heappop(pq)

        if visited[r][c]:
            continue
            
        # Calculate volume to reach this spill level
        if spill_level > water_level:
            volume_to_add = (spill_level - water_level) * len(pool_cells)
            if volume_to_add > 0:
                volume_increments.append(volume_to_add)
            total_volume += volume_to_add
            water_level = spill_level

        # If this is the target, we're done with calculation
        if (r, c) == target_node:
            break

        # Flood fill to find all newly connected cells at this water level
        q_flood = [(r, c)]
        visited[r][c] = True
        head = 0
        
        component_to_add = []
        while head < len(q_flood):
            curr_r, curr_c = q_flood[head]
            head += 1
            component_to_add.append((curr_r, curr_c))

            for dr, dc in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
                nr, nc = curr_r + dr, curr_c + dc
                if 0 <= nr < rows and 0 <= nc < cols and not visited[nr][nc]:
                    wall_height = max(depths[curr_r][curr_c], depths[nr][nc])
                    if wall_height <= water_level:
                        visited[nr][nc] = True
                        q_flood.append((nr, nc))

        # Add the new component to the pool and their neighbors to the pq
        for comp_r, comp_c in component_to_add:
             add_to_pool_and_pq(comp_r, comp_c)


    print("The total time is the sum of volumes needed to fill the pool to each subsequent spill level.")
    equation = " + ".join(map(str, volume_increments))
    print(f"Time = {equation}")
    print(f"Time = {total_volume} minutes")

solve_water_well()
<<<496>>>