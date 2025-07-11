import heapq

def solve_well_problem():
    """
    Calculates the time required for the water level in section 43 to begin to rise.
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
    source = (0, 0)
    target = (6, 6)

    # Step 2: Use Dijkstra's algorithm to find the min_barrier to each cell.
    # min_barrier[r][c] stores the minimum water level required to reach cell (r, c).
    min_barrier = [[float('inf')] * cols for _ in range(rows)]
    min_barrier[source[0]][source[1]] = depths[source[0]][source[1]]
    
    # Priority queue stores (max_depth_on_path, row, col)
    pq = [(depths[source[0]][source[1]], source[0], source[1])]

    while pq:
        d, r, c = heapq.heappop(pq)

        if d > min_barrier[r][c]:
            continue

        for dr, dc in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
            nr, nc = r + dr, c + dc
            
            if 0 <= nr < rows and 0 <= nc < cols:
                # The barrier to enter the neighbor is the max of the current path's barrier
                # and the neighbor's own depth.
                new_barrier = max(d, depths[nr][nc])
                
                if new_barrier < min_barrier[nr][nc]:
                    min_barrier[nr][nc] = new_barrier
                    heapq.heappush(pq, (new_barrier, nr, nc))

    # Step 3: Determine the final water level.
    spill_height = min_barrier[target[0]][target[1]]
    target_depth = depths[target[0]][target[1]]
    final_level = max(spill_height, target_depth)

    # Step 4: Calculate the total volume (time).
    total_volume = 0
    equation_terms = []

    for r in range(rows):
        for c in range(cols):
            # A cell is part of the final "lake" if it's reachable at the final_level.
            if min_barrier[r][c] <= final_level:
                volume_in_cell = final_level - depths[r][c]
                # We only sum non-zero volumes
                if volume_in_cell > 0:
                    total_volume += volume_in_cell
                    equation_terms.append(str(volume_in_cell))

    # Print the equation and the final answer as requested.
    equation_str = " + ".join(equation_terms)
    print(f"{equation_str} = {total_volume}")

solve_well_problem()
<<<1560>>>