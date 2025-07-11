import heapq

def solve_water_well_problem():
    """
    Solves the water well puzzle by finding the critical water level and calculating the total volume.
    """
    depths = [
        [1, 5, 27, 22, 28, 40, 14],
        [39, 13, 17, 30, 41, 12, 2],
        [32, 35, 24, 25, 19, 47, 34],
        [16, 33, 10, 42, 7, 44, 18],
        [3, 8, 45, 37, 4, 21, 20],
        [15, 46, 38, 6, 26, 48, 49],
        [9, 23, 31, 29, 11, 36, 43],
    ]
    rows, cols = 7, 7
    
    # Find start and target coordinates by their depths
    start_pos, target_pos = None, None
    for r in range(rows):
        for c in range(cols):
            if depths[r][c] == 1:
                start_pos = (r, c)
            if depths[r][c] == 43:
                target_pos = (r, c)

    # Part 1: Find the minimax path height using a modified Dijkstra's algorithm
    # Priority queue stores tuples of (max_depth_on_path, (row, col))
    pq = [(depths[start_pos[0]][start_pos[1]], start_pos)]
    min_max_depths = [[float('inf')] * cols for _ in range(rows)]
    min_max_depths[start_pos[0]][start_pos[1]] = depths[start_pos[0]][start_pos[1]]

    h_minimax = -1

    while pq:
        max_d, (r, c) = heapq.heappop(pq)

        if max_d > min_max_depths[r][c]:
            continue

        if (r, c) == target_pos:
            h_minimax = max_d
            break

        for dr, dc in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
            nr, nc = r + dr, c + dc

            if 0 <= nr < rows and 0 <= nc < cols:
                # The new max depth is the max of the path so far and the new cell's depth
                new_max_d = max(max_d, depths[nr][nc])
                
                if new_max_d < min_max_depths[nr][nc]:
                    min_max_depths[nr][nc] = new_max_d
                    heapq.heappush(pq, (new_max_d, (nr, nc)))
    
    # Part 2: Calculate the total volume
    # Find all cells reachable from start with depths <= h_minimax
    
    filled_sections = []
    total_volume = 0
    q = [start_pos]
    visited = {start_pos}

    while q:
        r, c = q.pop(0)
        filled_sections.append(depths[r][c])
        
        for dr, dc in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
            nr, nc = r + dr, c + dc
            if 0 <= nr < rows and 0 <= nc < cols and (nr, nc) not in visited:
                if depths[nr][nc] <= h_minimax:
                    visited.add((nr, nc))
                    q.append((nr, nc))

    # Build the equation string and calculate the sum
    equation_parts = []
    for depth in sorted(filled_sections):
        diff = h_minimax - depth
        total_volume += diff
        equation_parts.append(f"({h_minimax} - {depth})")
        
    final_equation = " + ".join(equation_parts)

    print(f"The minimum water level required for water to reach section 43 is {h_minimax}'.")
    print("This level is determined by the path from section 1 to 43 that has the lowest maximum depth.")
    print("\nThe total time required is the volume needed to fill all reachable sections up to this level.")
    print("The filled sections are those reachable from section 1 through paths with depths no greater than 43'.")
    print("\nThe volume calculation is:")
    print(f"{final_equation} = {total_volume}")
    print(f"\nAfter {total_volume} minutes will the water level on section 43 begin to rise.")


solve_water_well_problem()