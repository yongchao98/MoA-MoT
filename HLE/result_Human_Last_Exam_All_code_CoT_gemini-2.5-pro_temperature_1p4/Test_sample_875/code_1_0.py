import heapq

def solve_well_problem():
    """
    Solves the well-filling problem by determining the bottleneck water level,
    finding all connected cells at that level, and calculating the volume
    required to fill the resulting lake up to the target cell's bottom.
    """
    # The grid of depths provided in the problem
    grid_str = """
    1 5 27 22 28 40 14
    39 13 17 30 41 12 2
    32 35 24 25 19 47 34
    16 33 10 42 7 44 18
    3 8 45 37 4 21 20
    15 46 38 6 26 48 49
    9 23 31 29 11 36 43
    """
    grid = [list(map(int, row.split())) for row in grid_str.strip().split('\n')]
    rows, cols = len(grid), len(grid[0])

    # Create maps for easy lookup between depth and coordinates
    depth_to_coords = {}
    coords_to_depth = {}
    for r in range(rows):
        for c in range(cols):
            depth = grid[r][c]
            depth_to_coords[depth] = (r, c)
            coords_to_depth[(r, c)] = depth

    # --- Step 1: Find the bottleneck level to reach cell 43 ---
    # We use a modified Dijkstra's algorithm to find the minimax path.
    # The "distance" is the bottleneck level (maximum spill-over depth on the path).
    source_depth = 1
    target_depth = 43
    target_level = 43

    # dist[d] stores the minimum bottleneck level to reach depth d
    dist = {d: float('inf') for d in depth_to_coords}
    dist[source_depth] = source_depth

    # Priority queue stores (bottleneck_level, depth)
    pq = [(source_depth, source_depth)]
    
    bottleneck_level_to_reach_target = float('inf')

    while pq:
        max_level, current_d = heapq.heappop(pq)

        if max_level > dist[current_d]:
            continue

        if current_d == target_depth:
            # We found the best path to the target
            bottleneck_level_to_reach_target = max_level
            break
        
        r, c = depth_to_coords[current_d]
        
        # Check orthogonal neighbors
        for dr, dc in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
            nr, nc = r + dr, c + dc
            
            if 0 <= nr < rows and 0 <= nc < cols:
                neighbor_d = coords_to_depth[(nr, nc)]
                
                # Spill-over level is the depth of the shallower cell
                spill_level = min(current_d, neighbor_d)
                
                # The bottleneck for the path to the neighbor is the max of the
                # current path's bottleneck and the new spill-over level.
                new_bottleneck = max(max_level, spill_level)
                
                if new_bottleneck < dist[neighbor_d]:
                    dist[neighbor_d] = new_bottleneck
                    heapq.heappush(pq, (new_bottleneck, neighbor_d))

    # --- Step 2: Find all cells in the final lake ---
    # The lake consists of all cells reachable from the source using only
    # spill-over points at or below the determined bottleneck level.
    final_lake_cells_depths = set()
    queue = [source_depth]
    visited = {source_depth}

    while queue:
        current_d = queue.pop(0)
        final_lake_cells_depths.add(current_d)
        
        r, c = depth_to_coords[current_d]
        
        for dr, dc in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
            nr, nc = r + dr, c + dc
            
            if 0 <= nr < rows and 0 <= nc < cols:
                neighbor_d = coords_to_depth[(nr, nc)]
                spill_level = min(current_d, neighbor_d)
                
                if neighbor_d not in visited and spill_level <= bottleneck_level_to_reach_target:
                    visited.add(neighbor_d)
                    queue.append(neighbor_d)

    # --- Step 3: Calculate the total volume to fill the lake ---
    # The water level on section 43 rises when the lake's surface reaches its bottom (depth 43).
    # Total volume = Sum of (cell_depth - water_level_depth) for all cells deeper than the water level.
    volume_terms = []
    total_volume = 0
    for d in sorted(list(final_lake_cells_depths)): # sorted for consistent output
        if d > target_level:
            volume = d - target_level
            total_volume += volume
            volume_terms.append(f"({d} - {target_level})")

    # Final output
    print("The total volume is the sum of volumes for each section in the connected water body that is deeper than 43'.")
    print("The time in minutes for the water level on section 43 to begin to rise is:")
    if not volume_terms:
        print(f"0")
    else:
        equation = " + ".join(volume_terms)
        print(f"{equation} = {total_volume}")
    
    return total_volume

# Execute the solution
solve_well_problem()