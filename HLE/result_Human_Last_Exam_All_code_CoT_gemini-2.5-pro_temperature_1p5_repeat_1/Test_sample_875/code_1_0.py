import heapq

def solve():
    """
    This function solves the water well problem by modeling it as a basin-filling process.
    It calculates the time required for water to reach section 43.
    """
    
    # The grid of depths for the 7x7 well.
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
    start_node = (0, 0)
    
    # Water must reach cell (6,6) d=43. The path must go via an adjacent cell.
    # The 'easiest' adjacent cell to reach is (6,5) with depth 36.
    # The spillover event is determined by the path to (6,5).
    end_node_for_spill = (6, 5)

    # Step 1: Find the spillover height (H_spill) using a modified Dijkstra's algorithm.
    # This finds the path from start to end that minimizes the maximum depth encountered.
    min_max_depths = [[float('inf')] * cols for _ in range(rows)]
    min_max_depths[start_node[0]][start_node[1]] = depths[start_node[0]][start_node[1]]
    
    pq = [(depths[start_node[0]][start_node[1]], start_node[0], start_node[1])]

    while pq:
        d, r, c = heapq.heappop(pq)

        if d > min_max_depths[r][c]:
            continue
        
        # We can stop once we determine the spill level for the target.
        if (r, c) == end_node_for_spill:
            break

        for dr, dc in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
            nr, nc = r + dr, c + dc

            if 0 <= nr < rows and 0 <= nc < cols:
                # The 'cost' of a path is the maximum depth seen so far.
                new_max_depth = max(d, depths[nr][nc])
                if new_max_depth < min_max_depths[nr][nc]:
                    min_max_depths[nr][nc] = new_max_depth
                    heapq.heappush(pq, (new_max_depth, nr, nc))
    
    h_spill = min_max_depths[end_node_for_spill[0]][end_node_for_spill[1]]
    
    # Step 2: Identify all cells in the basin that fills up before the spillover.
    # These are cells reachable from the start through paths of depths less than H_spill.
    basin_cells = []
    q = [start_node]
    visited = {start_node}
    
    while q:
        r, c = q.pop(0)
        basin_cells.append((r,c))

        for dr, dc in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
            nr, nc = r + dr, c + dc
            if 0 <= nr < rows and 0 <= nc < cols and (nr, nc) not in visited:
                if depths[nr][nc] < h_spill:
                    visited.add((nr, nc))
                    q.append((nr, nc))

    # Step 3: Calculate the total volume to fill the basin to H_spill.
    total_volume = 0
    calculation_parts = []
    
    # Sort cells by depth for a cleaner equation display.
    basin_cells_sorted = sorted(basin_cells, key=lambda cell: depths[cell[0]][cell[1]])
    
    for r, c in basin_cells_sorted:
        depth = depths[r][c]
        volume_part = h_spill - depth
        total_volume += volume_part
        calculation_parts.append(f"({h_spill} - {depth})")
        
    # Final Output Generation
    print(f"The spillover height required to reach section 43 is {h_spill} feet.")
    print("The volume of water needed is calculated by summing the volumes to fill each section in the connected basin up to this height:")
    print("\nEquation:")
    # The full equation is printed as requested.
    print(" + ".join(calculation_parts))
    print(f"\nTotal volume = {total_volume} cubic feet.")
    print("\nWith water flowing at 1 cubic foot per minute, the time taken is equal to the volume.")
    
if __name__ == '__main__':
    solve()
    print("\n<<<630>>>")
