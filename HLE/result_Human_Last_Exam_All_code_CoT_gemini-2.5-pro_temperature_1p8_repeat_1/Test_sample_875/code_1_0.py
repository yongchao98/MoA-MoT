import heapq

def solve_water_well():
    """
    Calculates the time until water level rises in section 43.
    This is equivalent to the total volume of water poured.
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
    R, C = 7, 7
    start_pos = (0, 0)
    target_pos = (6, 6)

    # Water will enter the target from the adjacent cell that creates the lowest spill level.
    # Nbs of (6,6) are (6,5)[d=36] and (5,6)[d=49].
    # Spill level from (6,5) is max(36, 43) = 43.
    # Spill level from (5,6) is max(49, 43) = 49.
    # So water enters from (6,5). We need to flood (6,5) first.
    penultimate_target_pos = (6, 5)

    # Initial state
    h = depths[start_pos[0]][start_pos[1]]
    volume = 0
    basin = {start_pos}
    frontier = []  # A min-heap (priority queue)
    visited_frontier = {start_pos}
    volume_parts = []

    def add_neighbors_to_frontier(r, c):
        for dr, dc in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
            nr, nc = r + dr, c + dc
            if 0 <= nr < R and 0 <= nc < C and (nr, nc) not in visited_frontier:
                heapq.heappush(frontier, (depths[nr][nc], nr, nc))
                visited_frontier.add((nr, nc))

    add_neighbors_to_frontier(start_pos[0], start_pos[1])

    while frontier:
        d_spill, r_spill, c_spill = heapq.heappop(frontier)

        # Skip if this cell somehow got into the basin already
        if (r_spill, c_spill) in basin:
            continue
        
        # Check if we are about to flood the penultimate target
        if (r_spill, c_spill) == penultimate_target_pos:
            # Volume to reach the spill level of the penultimate target
            if h < d_spill:
                volume_added = (d_spill - h) * len(basin)
                volume_parts.append(volume_added)
                volume += volume_added
                h = d_spill
            
            # Now, calculate the final volume to spill into the target
            # The basin now includes the penultimate cell
            final_basin_size = len(basin) + 1
            final_h = depths[target_pos[0]][target_pos[1]]
            
            if h < final_h:
                final_volume_added = (final_h - h) * final_basin_size
                volume_parts.append(final_volume_added)
                volume += final_volume_added

            # We have the final answer, so we can stop.
            break

        # Calculate volume needed to fill the current basin to the spill level
        if h < d_spill:
            volume_added = (d_spill - h) * len(basin)
            volume_parts.append(volume_added)
            volume += volume_added
            h = d_spill

        # The new cell is now part of the basin
        basin.add((r_spill, c_spill))
        # Add its neighbors to the frontier
        add_neighbors_to_frontier(r_spill, c_spill)

    equation_str = ' + '.join(map(str, volume_parts))
    print(f"The total time is calculated by summing the volumes needed for each fill step:")
    print(f"{equation_str} = {volume}")
    return volume

final_answer = solve_water_well()
print(f"\nAfter {final_answer} minutes will the water level on section 43 begin to rise.")
print(f"<<<{final_answer}>>>")