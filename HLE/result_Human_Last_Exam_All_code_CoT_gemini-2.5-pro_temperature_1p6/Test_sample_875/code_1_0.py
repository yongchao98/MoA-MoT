import heapq

def solve_well_problem():
    """
    Calculates the time required for water to reach section 43.
    """
    # The grid of depths for the 7x7 well
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
    # The problem specifies the source above section 1, and the target is section 43
    source_pos = (0, 0)
    target_pos = (6, 6)

    # Step 1: Find the bottleneck water level `L` using a modified Dijkstra's algorithm.
    # This level is the minimum possible "maximum depth" on a path from source to target.
    min_max_costs = [[float('inf')] * cols for _ in range(rows)]
    r_s, c_s = source_pos
    min_max_costs[r_s][c_s] = depths[r_s][c_s]
    pq = [(depths[r_s][c_s], source_pos)] # (cost, (r, c))

    bottleneck_level = -1
    while pq:
        cost, (r, c) = heapq.heappop(pq)

        if cost > min_max_costs[r][c]:
            continue

        if (r, c) == target_pos:
            bottleneck_level = cost
            break

        for dr, dc in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
            nr, nc = r + dr, c + dc
            if 0 <= nr < rows and 0 <= nc < cols:
                new_cost = max(cost, depths[nr][nc])
                if new_cost < min_max_costs[nr][nc]:
                    min_max_costs[nr][nc] = new_cost
                    heapq.heappush(pq, (new_cost, (nr, nc)))

    # Step 2: Identify all cells in the basin at the bottleneck water level `L`.
    # These are all cells reachable from the source via cells with depth <= L.
    basin_cells = set()
    q = [source_pos]
    visited = {source_pos}
    if depths[source_pos[0]][source_pos[1]] <= bottleneck_level:
        basin_cells.add(source_pos)
        while q:
            r, c = q.pop(0)
            for dr, dc in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
                nr, nc = r + dr, c + dc
                if 0 <= nr < rows and 0 <= nc < cols and (nr, nc) not in visited:
                    if depths[nr][nc] <= bottleneck_level:
                        visited.add((nr, nc))
                        basin_cells.add((nr, nc))
                        q.append((nr, nc))

    # Step 3: Calculate the total volume of water needed.
    # Volume = sum(L - depth) for each cell in the basin.
    sum_of_depths_in_basin = sum(depths[r][c] for r, c in basin_cells)
    num_cells_in_basin = len(basin_cells)
    total_volume = num_cells_in_basin * bottleneck_level - sum_of_depths_in_basin

    # Print the explanation and final calculation
    print(f"The minimum water level (L) required to connect section 1 to section 43 is {bottleneck_level}'.")
    print(f"At this level, the water basin contains {num_cells_in_basin} sections.")
    print("The total volume is calculated by summing (L - depth) for every section in the basin.")
    print("\nThe final equation is:")
    print(f"Volume = (Number of sections in basin) * L - (Sum of depths in basin)")
    print(f"Volume = {num_cells_in_basin} * {bottleneck_level} - {sum_of_depths_in_basin}")
    print(f"Volume = {num_cells_in_basin * bottleneck_level} - {sum_of_depths_in_basin}")
    print(f"Volume = {total_volume}")
    print(f"\nThus, it will take {total_volume} minutes for the water level on section 43 to begin to rise.")

solve_well_problem()