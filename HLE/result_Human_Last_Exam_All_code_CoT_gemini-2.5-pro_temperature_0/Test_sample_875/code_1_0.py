import heapq

def solve_water_well_problem():
    """
    Calculates the time required for water to reach section 43 in the well.

    The problem is modeled as finding a path in a grid, where the cost of visiting
    a cell is its depth. We use a variation of Prim's algorithm to simulate the
    water filling the shallowest available adjacent sections first. A priority queue
    is used to efficiently select the next section to fill based on the minimum depth.

    The total time is the sum of the depths of all sections filled before the
    water begins to rise in the target section (depth 43).
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
    
    # Find the target node coordinates by its depth
    target_depth = 43
    target_node = None
    for r in range(rows):
        for c in range(cols):
            if depths[r][c] == target_depth:
                target_node = (r, c)
                break
        if target_node:
            break

    # The priority queue will store tuples of (depth, row, col)
    pq = []
    
    # Start the process at the source cell
    start_r, start_c = start_node
    heapq.heappush(pq, (depths[start_r][start_c], start_r, start_c))

    visited = [[False for _ in range(cols)] for _ in range(rows)]
    total_volume = 0
    path_of_depths = []

    while pq:
        # Get the next shallowest section from the frontier
        depth, r, c = heapq.heappop(pq)

        if visited[r][c]:
            continue

        visited[r][c] = True

        # If we have reached the target section, stop.
        # The accumulated volume is the time until this point.
        if (r, c) == target_node:
            break

        # Add the volume of the current section to the total
        total_volume += depth
        path_of_depths.append(depth)

        # Add its unvisited neighbors to the frontier
        for dr, dc in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
            nr, nc = r + dr, c + dc

            if 0 <= nr < rows and 0 <= nc < cols and not visited[nr][nc]:
                heapq.heappush(pq, (depths[nr][nc], nr, nc))

    # Format and print the final equation and result
    equation = " + ".join(map(str, path_of_depths))
    print("The water fills sections with the following depths in order:")
    print(f"{equation} = {total_volume}")
    print(f"\nAfter {total_volume} minutes, the water level on section {target_depth} will begin to rise.")

solve_water_well_problem()