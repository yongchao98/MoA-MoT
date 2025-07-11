import heapq

def solve_water_well_problem():
    """
    Solves the water well filling problem by:
    1. Finding the minimum water level required to reach the target cell using a Dijkstra-like algorithm.
    2. Identifying all cells in the basin that fills up to this level.
    3. Calculating the total volume of water required by summing the volume in each basin cell.
    4. Printing the step-by-step calculation and the final answer.
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
    start_pos = (0, 0)
    end_pos = (6, 6)

    # Step 1: Find the final water level using a Dijkstra-like algorithm
    # to find the path with the minimum bottleneck (minimum max depth).
    min_levels = [[float('inf')] * cols for _ in range(rows)]
    pq = [(depths[start_pos[0]][start_pos[1]], start_pos[0], start_pos[1])]
    min_levels[start_pos[0]][start_pos[1]] = depths[start_pos[0]][start_pos[1]]

    final_level = -1

    while pq:
        level, r, c = heapq.heappop(pq)

        if level > min_levels[r][c]:
            continue

        if (r, c) == end_pos:
            final_level = level
            break

        for dr, dc in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
            nr, nc = r + dr, c + dc

            if 0 <= nr < rows and 0 <= nc < cols:
                # The new water level is the max of the current level and the depth of the new cell
                new_level = max(level, depths[nr][nc])
                if new_level < min_levels[nr][nc]:
                    min_levels[nr][nc] = new_level
                    heapq.heappush(pq, (new_level, nr, nc))
    
    # Step 2 & 3: Identify basin cells and calculate total volume
    basin_cells_depths = []
    for r in range(rows):
        for c in range(cols):
            # A cell is in the basin if its depth is less than or equal to the final spill-over level.
            # We assume all such cells are connected, which is true for this grid.
            if depths[r][c] <= final_level:
                basin_cells_depths.append(depths[r][c])

    total_volume = 0
    calculation_parts = []
    
    # Sort depths for a more readable equation
    basin_cells_depths.sort()

    for depth in basin_cells_depths:
        volume_in_cell = final_level - depth
        total_volume += volume_in_cell
        calculation_parts.append(f"({final_level} - {depth})")

    # Step 4: Print the formatted output
    print(f"To reach section 43, the water level must rise to a depth of {final_level}'.")
    print("The basin consists of all connected sections with a depth less than or equal to this level.")
    print("The total time required is the total volume of water needed to fill this basin up to {final_level}'.")
    print("\nThe calculation for the total volume is:")
    
    equation_str = " + ".join(calculation_parts)
    print(f"{equation_str} = {total_volume}")

    print(f"\nTherefore, the water level on section 43 will begin to rise after {total_volume} minutes.")


solve_water_well_problem()
<<<903>>>