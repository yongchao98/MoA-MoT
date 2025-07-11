import collections

def solve_maze():
    """
    Solves the maze problem by finding the least dangerous path from '@' to 'g'.
    The "least dangerous" path is interpreted as the shortest path that does not
    go through any dangerous tiles (the dragon 'D').
    """
    # The map is parsed from the multi-line string into a 2D grid.
    # Note: The provided LaTeX-like format is converted to a character grid.
    # I'm interpreting the grid based on character positions from the raw string.
    
    map_str = r"""/ &  & - & - & - & - & - & - &  &  &  &  &  &  &  &  &  &  &  &  &  &  &  &  &  &  &  &  &  &  \\
/ &  & | & . & . & . & . & | &  &  &  &  &  &  & \# & \# & \# & \# & \# & \# & \# & \# & \# & \# & \# & \# &  &  &  &  \\
/ &  & | & . & . & . & . & | &  &  &  &  &  &  & \# &  &  &  &  &  &  &  &  &  &  & \# &  &  &  &  \\
/ &  & | & . & g & . & . & + & \# & \# & \# & \# & \# & \# & \# & \# &  &  &  &  &  &  &  &  &  & @ &  &  &  &  \\
/ &  & | & . & . & . & . & | &  &  &  &  &  &  &  & \# &  &  &  &  &  &  & - & - & - & + & - & - & - &  \\
/ &  & - & - & - & - & - & - &  &  &  &  &  &  &  & \# &  &  &  &  &  &  & | & . & . & . & . & . & | \\
/ &  &  &  &  &  &  &  &  &  &  &  &  &  &  & \# &  &  &  &  &  &  & | & . & ! & . & . & . & | \\
/ &  &  &  &  &  &  &  &  &  &  &  &  &  &  & \# &  &  &  &  &  &  & | & . & . & . & . & . & | \\
/ &  &  &  &  &  &  &  &  &  &  &  &  &  &  & \# &  &  &  &  &  &  & | & . & . & . & . & . & | \\
/ &  &  &  & - & - & - & - &  &  &  &  &  &  &  & \# &  &  &  &  &  &  & | & . & . & . & . & . & | \\
/ &  &  &  & | & . & . & | &  &  &  &  &  &  &  & \# & \# & \# & \# & \# & \# & \# & + & . & . & D & . & . & | \\
/ &  &  &  & | & < & . & + & \# & \# & \# &  &  &  &  & \# &  &  &  &  &  &  & | & . & . & . & . & . & | \\
/ &  &  &  & - & - & - & - &  &  & \# &  &  &  &  & \# &  &  &  &  &  &  & | & . & ? & . & . & . & | \\
/ &  &  &  &  &  &  &  &  &  & \# & \# & \# & \# & \# & \# &  &  &  &  &  &  & - & - & - & - & - & - & - \\"""

    grid = [line.strip() for line in map_str.strip().split('\\')]

    # According to the user's key:
    # Walkable: '#', '.', '+', '!', '<', '?', 'g', '@', ' '
    # Walls/Dangers: '-', '/', '|', 'D'
    walls = {'-', '/', '|', 'D'}

    rows, cols = len(grid), max(len(row) for row in grid)
    
    start_pos, goal_pos = None, None
    for r in range(rows):
        for c in range(len(grid[r])):
            if grid[r][c] == '@':
                start_pos = (r, c)
            elif grid[r][c] == 'g':
                goal_pos = (r, c)

    if not start_pos or not goal_pos:
        print("Start or goal not found in map.")
        return

    # BFS implementation
    queue = collections.deque([(start_pos, [])])
    visited = {start_pos}
    
    path_to_goal = None

    while queue:
        (r, c), path = queue.popleft()

        if (r, c) == goal_pos:
            path_to_goal = path
            break

        # Explore neighbors: Up, Down, Left, Right
        for dr, dc, move in [(0, 1, 'R'), (0, -1, 'L'), (1, 0, 'D'), (-1, 0, 'U')]:
            nr, nc = r + dr, c + dc

            if 0 <= nr < rows and 0 <= nc < len(grid[nr]) and \
               grid[nr][nc] not in walls and (nr, nc) not in visited:
                
                visited.add((nr, nc))
                new_path = path + [move]
                queue.append(((nr, nc), new_path))
    
    if path_to_goal:
        # Compress the path: RRRLL -> RL
        if not path_to_goal:
            compressed_path = ""
        else:
            compressed_path_list = [path_to_goal[0]]
            for move in path_to_goal[1:]:
                if move != compressed_path_list[-1]:
                    compressed_path_list.append(move)
            compressed_path = "".join(compressed_path_list)
        print(compressed_path)
    else:
        print("No path found.")

solve_maze()