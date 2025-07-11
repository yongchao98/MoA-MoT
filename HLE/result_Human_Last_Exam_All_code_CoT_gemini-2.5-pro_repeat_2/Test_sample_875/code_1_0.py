import heapq

def solve():
    """
    Simulates water filling the well to find when the level at section 43 starts to rise.
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
    
    target_cell = (6, 6)
    rows, cols = 7, 7

    lake_cells = set()
    # Priority queue stores (depth, r, c) for cells on the frontier
    # We start with the source cell (0, 0)
    frontier_pq = [(depths[0][0], 0, 0)]
    visited = {(0, 0)}
    
    total_time = 0
    current_level = 0
    
    time_components = []

    while frontier_pq:
        # Get the next spill event from the frontier
        spill_level, r, c = heapq.heappop(frontier_pq)

        # This check handles cases where a cell is added to the frontier but later gets flooded
        # by a higher water level before becoming a spill point itself.
        if (r, c) in lake_cells:
            continue

        # Calculate time to fill the lake to the new spill_level
        if lake_cells:
            time_added = (spill_level - current_level) * len(lake_cells)
            if time_added > 0:
                time_components.append(f"({spill_level} - {current_level}) * {len(lake_cells)}")
                total_time += time_added
        else:
            # Special case for the very first cell
            time_added = spill_level
            time_components.append(str(time_added))
            total_time += time_added
            
        current_level = spill_level

        # Expand the lake using a flood-fill (BFS) from the spill gateway
        # This finds all newly connected cells at or below the current water level
        q = [(r, c)]
        visited.add((r, c))
        
        newly_added_to_lake = []
        
        while q:
            curr_r, curr_c = q.pop(0)
            
            # Check if target is reached
            if (curr_r, curr_c) == target_cell:
                print("Calculation steps:")
                print(" + ".join(time_components))
                print(f"\nTotal time = {' + '.join(str(eval(comp)) for comp in time_components)} = {total_time}")
                print(f"\nThe water level on section 43 will begin to rise after {total_time} minutes.")
                return

            lake_cells.add((curr_r, curr_c))
            newly_added_to_lake.append((curr_r, curr_c))

            # Check neighbors for flooding or adding to the new frontier
            for dr, dc in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
                nr, nc = curr_r + dr, curr_c + dc
                
                if 0 <= nr < rows and 0 <= nc < cols and (nr, nc) not in visited:
                    neighbor_depth = depths[nr][nc]
                    if neighbor_depth <= current_level:
                        q.append((nr, nc))
                    else:
                        heapq.heappush(frontier_pq, (neighbor_depth, nr, nc))
                    visited.add((nr, nc))
        
        # After flooding, add neighbors of all newly added lake cells to the frontier
        for lake_r, lake_c in newly_added_to_lake:
            for dr, dc in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
                nr, nc = lake_r + dr, lake_c + dc
                if 0 <= nr < rows and 0 <= nc < cols and (nr, nc) not in visited:
                    heapq.heappush(frontier_pq, (depths[nr][nc], nr, nc))
                    visited.add((nr, nc))

solve()
<<<333>>>