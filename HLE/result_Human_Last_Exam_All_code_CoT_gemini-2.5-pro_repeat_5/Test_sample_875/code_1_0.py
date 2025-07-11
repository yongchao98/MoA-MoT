def solve_well_problem():
    """
    Calculates the time required for the water level on section 43 to begin to rise.
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

    target_depth = 43
    rows, cols = len(depths), len(depths[0])

    # Find the coordinates of the deepest section
    max_depth = 0
    start_pos = None
    for r in range(rows):
        for c in range(cols):
            if depths[r][c] > max_depth:
                max_depth = depths[r][c]
                start_pos = (r, c)

    # Find the connected component of sections deeper than the target_depth,
    # starting from the deepest section. This identifies the relevant basin.
    # We use Breadth-First Search (BFS).
    basin_sections = []
    if start_pos and depths[start_pos[0]][start_pos[1]] > target_depth:
        q = [start_pos]
        visited = {start_pos}
        basin_sections.append(start_pos)
        
        while q:
            r, c = q.pop(0)
            
            for dr, dc in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
                nr, nc = r + dr, c + dc
                
                if 0 <= nr < rows and 0 <= nc < cols and (nr, nc) not in visited:
                    if depths[nr][nc] > target_depth:
                        visited.add((nr, nc))
                        q.append((nr, nc))
                        basin_sections.append((nr, nc))

    # Calculate the volume needed to fill this basin up to the target_depth level
    total_volume = 0
    equation_parts = []
    
    # Sort for consistent output order
    basin_sections.sort(key=lambda p: depths[p[0]][p[1]])

    for r, c in basin_sections:
        depth = depths[r][c]
        volume_contribution = depth - target_depth
        total_volume += volume_contribution
        equation_parts.append(f"({depth} - {target_depth})")

    print("The water level on section 43 begins to rise when the basin it belongs to fills up to its bottom.")
    print(f"The target water level corresponds to the bottom of section 43, which is at a depth of {target_depth}'.")
    print("We must calculate the volume of all sections in the deepest basin that are deeper than 43'.")
    print("The relevant sections in the deepest basin are those with depths:", ", ".join(str(depths[r][c]) for r, c in basin_sections))
    print("The calculation for the total volume is the sum of (depth - 43) for each of these sections:")
    print(f"Volume = {' + '.join(equation_parts)}")
    
    calculation_result = " + ".join(str(depths[r][c] - target_depth) for r,c in basin_sections)
    print(f"Volume = {calculation_result}")
    
    print(f"Total Volume = {total_volume} cubic feet.")
    print("Since the water is poured at a rate of 1 cubic foot per minute, the time is equal to the volume.")
    print(f"Time = {total_volume} minutes.")

solve_well_problem()