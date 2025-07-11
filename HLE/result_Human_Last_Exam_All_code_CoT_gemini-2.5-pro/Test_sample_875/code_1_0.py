import heapq

def solve():
    """
    Solves the water well problem by determining the volume of water required
    for the level to start rising in the target cell.
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
    end_node = (6, 6)

    # Step 1: Find the critical depth using a modified Dijkstra's algorithm.
    # This finds the path from start to end that minimizes the maximum depth encountered.
    min_max_depths = [[float('inf')] * cols for _ in range(rows)]
    # The priority queue stores tuples of (max_depth_on_path, row, col).
    pq = [(depths[start_node[0]][start_node[1]], start_node[0], start_node[1])]
    min_max_depths[start_node[0]][start_node[1]] = depths[start_node[0]][start_node[1]]

    while pq:
        max_depth, r, c = heapq.heappop(pq)

        if max_depth > min_max_depths[r][c]:
            continue

        # If we've reached the end, this is the best path's bottleneck.
        if (r, c) == end_node:
            break

        # Explore neighbors
        for dr, dc in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
            nr, nc = r + dr, c + dc

            if 0 <= nr < rows and 0 <= nc < cols:
                # The new path's bottleneck is the max of the current path's bottleneck
                # and the depth of the new cell.
                new_max_depth = max(max_depth, depths[nr][nc])
                
                if new_max_depth < min_max_depths[nr][nc]:
                    min_max_depths[nr][nc] = new_max_depth
                    heapq.heappush(pq, (new_max_depth, nr, nc))

    critical_depth = min_max_depths[end_node[0]][end_node[1]]

    # Step 2: Identify all cells in the basin using BFS.
    # The basin includes all cells reachable from the start without exceeding the critical_depth.
    basin_cells = set()
    q = [start_node]
    visited = {start_node}

    while q:
        r, c = q.pop(0)
        
        # A cell is part of the basin if it's connected and its depth is not a "wall".
        if depths[r][c] <= critical_depth:
            basin_cells.add((r, c))
            for dr, dc in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
                nr, nc = r + dr, c + dc
                if 0 <= nr < rows and 0 <= nc < cols and (nr, nc) not in visited:
                    visited.add((nr, nc))
                    q.append((nr, nc))

    # Step 3: Calculate the total volume and format the output.
    total_volume = 0
    # Sort cells for a consistent and readable output order.
    sorted_basin_cells = sorted(list(basin_cells))

    equation_str_parts = []
    for r, c in sorted_basin_cells:
        # The volume in each cell is the difference between the water level and the cell's bottom.
        volume_in_cell = critical_depth - depths[r][c]
        total_volume += volume_in_cell
        equation_str_parts.append(f"({critical_depth} - {depths[r][c]})")

    print(f"The minimum water level required to reach section 43 is {critical_depth} feet.")
    print(f"The water basin at this level contains {len(basin_cells)} cells.")
    print("\nThe total volume is the sum of the water required to fill each cell in the basin up to this level.")
    print(f"Volume = {' + '.join(equation_str_parts)}")
    print(f"\nTotal Volume = {total_volume} cubic feet.")
    print(f"Therefore, the time for the water level on section 43 to begin to rise is {total_volume} minutes.")

solve()