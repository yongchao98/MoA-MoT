import collections

def solve_maze():
    """
    Solves the maze problem by parsing the map, running BFS to find the
    shortest path from '@' to 'g' while avoiding 'D', and then formatting
    the output path.
    """
    raw_map_data = r"""
    / &  & - & - & - & - & - & - &  &  &  &  &  &  &  &  &  &  &  &  &  &  &  &  &  &  &  &  &  &  \\
    / &  & | & . & . & . & . & | &  &  &  &  &  &  & \# & \# & \# & \# & \# & \# & \# & \# & \# & \# & \# & \# &  &  &  &  \\
    / &  & | & . & . & . & . & | &  &  &  &  &  &  & \# &  &  &  &  &  &  &  &  &  &  & \# &  &  &  &  \\
    / &  & | & . & g & . & . & + & \# & \# & \# & \# & \# & \# & \# & \# &  &  &  &  &  &  &  &  &  & @ &  &  &  &  \\
    / &  & | & . & . & . & . & | &  &  &  &  &  &  &  & \# &  &  &  &  &  &  & - & - & - & + & - & - & - &  \\
    / &  & - & - & - & - & - & - &  &  &  &  &  &  &  & \# &  &  &  &  &  &  & | & . & . & . & . & . & | \\
    / &  &  &  &  &  &  &  &  &  &  &  &  &  &  & \# &  &  &  &  &  &  & | & . & ! & . & . & . & | \\
    / &  &  &  &  &  &  &  &  &  &  &  &  &  &  & \# &  &  &  &  &  &  & | & . & . & . & . & . & | \\
    / &  &  &  &  &  &  &  &  &  &  &  &  &  &  & \# &  &  &  &  &  &  & | & . & . & . & . & . & | \\
    / &  &  &  & - & - & - & - &  &  &  &  &  &  &  & \# &  &  &  &  &  &  & | & . & . & . & . & | \\
    / &  &  &  & | & . & . & | &  &  &  &  &  &  &  & \# & \# & \# & \# & \# & \# & \# & + & . & . & D & . & . & | \\
    / &  &  &  & | & < & . & + & \# & \# & \# &  &  &  &  & \# &  &  &  &  &  &  & | & . & . & . & . & . & | \\
    / &  &  &  & - & - & - & - &  &  & \# & \# & \# & \# & \# & \# &  &  &  &  &  &  & | & . & ? & . & . & . & | \\
    / &  &  &  &  &  &  &  &  &  & \# & \# & \# & \# & \# & \# &  &  &  &  &  &  & - & - & - & - & - & - & - \\
    """

    # 1. Parse the map from the raw string into a 2D grid
    lines = raw_map_data.strip().split('\\\\')
    grid_from_split = []
    for line in lines:
        line = line.strip()
        if not line: continue
        grid_from_split.append([c.strip() for c in line.split('&')])

    grid = []
    max_len = 0
    for row_list in grid_from_split:
        row = [item[0] if item else ' ' for item in row_list]
        grid.append(row)
        if len(row) > max_len:
            max_len = len(row)

    for row in grid:
        while len(row) < max_len:
            row.append(' ')

    # 2. Find start, goal, and define walls
    start_pos, goal_pos = None, None
    rows, cols = len(grid), max_len
    for r in range(rows):
        for c in range(cols):
            if grid[r][c] == '@':
                start_pos = (r, c)
            elif grid[r][c] == 'g':
                goal_pos = (r, c)
    
    # Define impassable tiles for the pathfinding algorithm.
    # The Dragon 'D' is treated as a wall to find the "least dangerous" path.
    # Empty space ' ' is treated as a wall to be safe.
    walls = {'#', '/', ' ', 'D'}

    # 3. BFS to find the shortest path
    queue = collections.deque([(start_pos, [])])
    visited = {start_pos}
    final_path = None

    while queue:
        (r, c), path = queue.popleft()

        if (r, c) == goal_pos:
            final_path = path
            break

        # Explore neighbors in Up, Down, Left, Right order
        moves = [(-1, 0, 'U'), (1, 0, 'D'), (0, -1, 'L'), (0, 1, 'R')]
        for dr, dc, move_char in moves:
            nr, nc = r + dr, c + dc

            if not (0 <= nr < rows and 0 <= nc < cols) or \
               (nr, nc) in visited or \
               grid[nr][nc] in walls:
                continue
            
            visited.add((nr, nc))
            new_path = path + [move_char]
            queue.append(((nr, nc), new_path))

    # 4. Condense the path format
    if not final_path:
        print("No path found.")
        return

    condensed_path = []
    if final_path:
        last_move = None
        for move in final_path:
            if move != last_move:
                condensed_path.append(move)
                last_move = move

    print("".join(condensed_path))

solve_maze()