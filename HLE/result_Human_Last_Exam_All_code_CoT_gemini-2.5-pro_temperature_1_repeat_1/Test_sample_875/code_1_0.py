import heapq

def solve():
    """
    Solves the water well problem by simulating the filling process.
    The simulation is based on a shortest path algorithm (Dijkstra) on a grid,
    which determines the sequence of spillover events.
    """
    # The depths of the 7x7 well sections.
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
    target_val = 43

    # Find the position of the target section.
    target_pos = None
    for r in range(rows):
        for c in range(cols):
            if depths[r][c] == target_val:
                target_pos = (r, c)
                break
    
    target_depth = depths[target_pos[0]][target_pos[1]]

    # Step 1: Use a Dijkstra-like algorithm to find the minimum water level
    # required to flood each cell, starting from the source at (0,0).
    # The 'distance' to a cell is the max depth on the path of least resistance.
    min_level = [[float('inf')] * cols for _ in range(rows)]
    pq = [(depths[start_pos[0]][start_pos[1]], start_pos[0], start_pos[1])]
    min_level[start_pos[0]][start_pos[1]] = depths[start_pos[0]][start_pos[1]]

    while pq:
        level, r, c = heapq.heappop(pq)

        if level > min_level[r][c]:
            continue

        for dr, dc in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
            nr, nc = r + dr, c + dc
            if 0 <= nr < rows and 0 <= nc < cols:
                # To flow to a neighbor, the water level must rise to the higher of
                # the current cell's spill-level and the neighbor's floor depth.
                new_level = max(level, depths[nr][nc])
                if new_level < min_level[nr][nc]:
                    min_level[nr][nc] = new_level
                    heapq.heappush(pq, (new_level, nr, nc))

    # Step 2: Calculate the total volume (time) by simulating the filling process
    # using the computed spillover levels.
    
    # Get all unique spillover levels from the grid in sorted order.
    # These are the critical water levels where the lake expands.
    all_levels = sorted(list(set(l for row in min_level for l in row if l != float('inf'))))

    total_volume = 0
    current_level = all_levels[0]
    equation_parts = []

    for i in range(1, len(all_levels)):
        # Determine the size of the lake at the current water level.
        # The lake consists of all cells that can be flooded at this level.
        lake_size = 0
        for r in range(rows):
            for c in range(cols):
                if min_level[r][c] <= current_level:
                    lake_size += 1
        
        next_level = all_levels[i]

        # Check if the target cell's floor is reached in this filling step.
        if current_level < target_depth <= next_level:
            volume_to_add = (target_depth - current_level) * lake_size
            total_volume += volume_to_add
            equation_parts.append(f"({target_depth} - {current_level}) * {lake_size}")
            # The answer is found. Print the full equation and the result.
            print(" + ".join(equation_parts) + f" = {int(total_volume)}")
            return int(total_volume)

        # If target not reached yet, fill up to the next spillover level.
        volume_to_add = (next_level - current_level) * lake_size
        total_volume += volume_to_add
        equation_parts.append(f"({next_level} - {current_level}) * {lake_size}")
        
        current_level = next_level

    return -1 # Should not be reached if target is in the grid.

final_answer = solve()
print(f"<<<{final_answer}>>>")