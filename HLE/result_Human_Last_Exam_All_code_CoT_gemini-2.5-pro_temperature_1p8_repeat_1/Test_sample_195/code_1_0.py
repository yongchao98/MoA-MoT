import collections

def solve_dungeon_path():
    """
    Finds the least dangerous path from the adventurer '@' to the gold 'g'
    by parsing the map, running a Breadth-First Search to find the shortest
    path while avoiding the dragon 'D', and then condensing the path string.
    """
    map_str = r"""
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
    / &  &  &  & - & - & - & - &  &  & \# &  &  &  &  & \# &  &  &  &  &  &  & | & . & ? & . & . & . & | \\
    / &  &  &  &  &  &  &  &  &  & \# & \# & \# & \# & \# & \# &  &  &  &  &  &  & - & - & - & - & - & - & - \\
    """

    # 1. Parse the map from bmatrix format into a 2D grid
    rows_str = map_str.strip().replace('\n', '').split(r'\\')
    grid = []
    for row_str in rows_str:
        if row_str.strip():
            cols = [c.strip() for c in row_str.split('&')]
            grid.append(cols)

    # 2. Find start, goal, and define obstacles
    start_pos, goal_pos = None, None
    for r, row in enumerate(grid):
        for c, char in enumerate(row):
            if char == '@':
                start_pos = (r, c)
            elif char == 'g':
                goal_pos = (r, c)
    
    if not start_pos or not goal_pos:
        print("Error: Start ('@') or goal ('g') not found on the map.")
        return

    # 'D' (Dragon) is dangerous, so it's a wall.
    walls = {'#', '-', '|', '/', 'D'}

    # 3. Use Breadth-First Search (BFS) to find the shortest path
    queue = collections.deque([(start_pos, "")])  # Each item: (position, path_string)
    visited = {start_pos}
    
    grid_rows, grid_cols = len(grid), len(grid[0])
    moves = {'U': (-1, 0), 'D': (1, 0), 'L': (0, -1), 'R': (0, 1)}

    path_to_goal = None
    
    while queue:
        (r, c), path = queue.popleft()
        
        if (r, c) == goal_pos:
            path_to_goal = path
            break

        for move_char, (dr, dc) in moves.items():
            nr, nc = r + dr, c + dc
            
            # Check if neighbor is within grid bounds and is a valid move
            if 0 <= nr < grid_rows and 0 <= nc < grid_cols:
                if (nr, nc) not in visited and grid[nr][nc] not in walls:
                    visited.add((nr, nc))
                    queue.append(((nr, nc), path + move_char))

    # 4. Condense the path string and print the result
    if path_to_goal:
        condensed_path = ""
        if len(path_to_goal) > 0:
            condensed_path += path_to_goal[0]
            for i in range(1, len(path_to_goal)):
                if path_to_goal[i] != path_to_goal[i-1]:
                    condensed_path += path_to_goal[i]
        
        print(condensed_path)
    else:
        print("No path to the gold could be found.")

solve_dungeon_path()
>>> DLDULUL