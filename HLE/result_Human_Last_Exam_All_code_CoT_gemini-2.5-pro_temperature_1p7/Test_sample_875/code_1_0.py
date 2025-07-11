import heapq

def solve_water_well_problem():
    """
    Solves the water well filling problem by:
    1. Modeling the well as a grid and the flow as filling interconnected basins.
    2. Using a Dijkstra-like algorithm to find the minimum spill level required to connect the source
       to the section that will spill into the target.
    3. Calculating the total volume of water needed to reach this spill level.
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

    N = 7
    start_node = (0, 0)
    # Water must flow from a shallower neighbor into the pit at (6,6), depth 43.
    # Neighbors of (6,6) are (5,6) [depth 49] and (6,5) [depth 36].
    # Water must come from (6,5).
    target_predecessor = (6, 5)

    # Step 1: Find the minimum spill level using a Dijkstra-like algorithm.
    # The 'distance' or 'cost' of a path is the maximum depth encountered.
    max_depths = [[float('inf')] * N for _ in range(N)]
    # The priority queue stores tuples of (current_max_depth, row, col).
    pq = [(depths[start_node[0]][start_node[1]], start_node[0], start_node[1])]
    max_depths[start_node[0]][start_node[1]] = depths[start_node[0]][start_node[1]]

    while pq:
        d, r, c = heapq.heappop(pq)

        if d > max_depths[r][c]:
            continue

        # Explore orthogonal neighbors.
        for dr, dc in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
            nr, nc = r + dr, c + dc

            if 0 <= nr < N and 0 <= nc < N:
                # The max depth of the path to the neighbor is the greater of the
                # current path's max depth and the neighbor's own depth.
                new_max_depth = max(d, depths[nr][nc])

                # If we found a "cheaper" path to the neighbor, update it.
                if new_max_depth < max_depths[nr][nc]:
                    max_depths[nr][nc] = new_max_depth
                    heapq.heappush(pq, (new_max_depth, nr, nc))

    # Step 2: The spill level is the minimum max_depth required to reach the predecessor.
    spill_level = max_depths[target_predecessor[0]][target_predecessor[1]]

    # Step 3: Calculate the volume required to fill the basin.
    # The basin consists of all cells reachable with a max_depth <= spill_level.
    basin_cells_info = []
    total_volume = 0
    for r in range(N):
        for c in range(N):
            if max_depths[r][c] <= spill_level:
                volume_for_cell = spill_level - depths[r][c]
                total_volume += volume_for_cell
                basin_cells_info.append({"depth": depths[r][c], "volume": volume_for_cell})

    # Step 4: Output the results as requested.
    basin_cells_info.sort(key=lambda x: x['depth'])
    
    print(f"The water level must rise to a height equivalent to depth {spill_level}.")
    print("This happens when the basin connecting the source to the section at (6,5) is full.")
    print("\nThe total time is the volume of water needed to fill all cells in this basin up to that level.")
    print("The final equation for the total volume is the sum of volumes for each cell in the basin:")
    
    # "output each number in the final equation"
    sum_parts = [f"({spill_level} - {cell['depth']})" for cell in basin_cells_info]
    print(f"Time = {' + '.join(sum_parts)}")
    
    # Also show the numeric values of each term
    volume_parts = [str(cell['volume']) for cell in basin_cells_info]
    print(f"Time = {' + '.join(volume_parts)}")
    
    print(f"Total Time = {total_volume} minutes")

solve_water_well_problem()