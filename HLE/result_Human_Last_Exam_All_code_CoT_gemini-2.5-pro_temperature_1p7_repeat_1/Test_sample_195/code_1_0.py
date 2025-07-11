import collections

def solve_dungeon_path():
    """
    This function solves the dungeon pathfinding problem.
    It parses the map, finds the start and end points,
    and uses Breadth-First Search (BFS) to find the shortest
    and safest path. The final path is simplified and printed.
    """
    
    # The map is represented as a raw string from the problem.
    # It will be parsed into a 2D grid.
    matrix_str = r"""
     &  & - & - & - & - & - & - &  &  &  &  &  &  &  &  &  &  &  &  &  &  &  &  &  &  &  &  &  \\
     &  & | & . & . & . & . & | &  &  &  &  &  &  & \# & \# & \# & \# & \# & \# & \# & \# & \# & \# & \# & \# &  &  &  \\
     &  & | & . & . & . & . & | &  &  &  &  &  &  & \# &  &  &  &  &  &  &  &  &  &  & \# &  &  &  \\
     &  & | & . & g & . & . & + & \# & \# & \# & \# & \# & \# & \# & \# &  &  &  &  &  &  &  &  &  & @ &  &  &  \\
     &  & | & . & . & . & . & | &  &  &  &  &  &  &  & \# &  &  &  &  &  &  & - & - & - & + & - & - & - \\
     &  & - & - & - & - & - & - &  &  &  &  &  &  &  & \# &  &  &  &  &  &  & | & . & . & . & . & . & | \\
     &  &  &  &  &  &  &  &  &  &  &  &  &  &  & \# &  &  &  &  &  &  & | & . & ! & . & . & . & | \\
     &  &  &  &  &  &  &  &  &  &  &  &  &  &  & \# &  &  &  &  &  &  & | & . & . & . & . & . & | \\
     &  &  &  &  &  &  &  &  &  &  &  &  &  &  & \# &  &  &  &  &  &  & | & . & . & . & . & . & | \\
     &  &  &  & - & - & - & - &  &  &  &  &  &  &  & \# &  &  &  &  &  &  & | & . & . & . & . & . & | \\
     &  &  &  & | & . & . & | &  &  &  &  &  &  &  & \# & \# & \# & \# & \# & \# & \# & + & . & . & D & . & . & | \\
     &  &  &  & | & < & . & + & \# & \# & \# &  &  &  &  & \# &  &  &  &  &  &  & | & . & . & . & . & . & | \\
     &  &  &  & - & - & - & - &  &  & \# &  &  &  &  & \# &  &  &  &  &  &  & | & . & ? & . & . & . & | \\
     &  &  &  &  &  &  &  &  &  & \# & \# & \# & \# & \# & \# &  &  &  &  &  &  & - & - & - & - & - & - & - \\
    """
    
    # Parse the string into a list of lists (grid)
    lines = matrix_str.strip().split(r"\\")
    grid = []
    for line in lines:
        line = line.strip()
        if not line:
            continue
        cells = [c.strip() if c.strip() else ' ' for c in line.split('&')]
        grid.append(cells)

    # Find the positions of the Adventurer (@) and the Gold (g)
    start_pos, goal_pos = None, None
    for r, row in enumerate(grid):
        for c, char in enumerate(row):
            if char == '@':
                start_pos = (r, c)
            elif char == 'g':
                goal_pos = (r, c)

    if not start_pos or not goal_pos:
        print("Error: Start ('@') or goal ('g') not found in the map.")
        return

    # Define tiles that are considered walkable. The Dragon 'D' is avoided for safety.
    walkable_chars = {'.', 'g', '+', '!', '@', '<', '?'}
    rows, cols = len(grid), len(grid[0])
    
    # Initialize BFS queue with starting position and an empty path
    queue = collections.deque([(start_pos, [])])  # State: ((row, col), path_list)
    visited = {start_pos}
    
    # Define the four cardinal directions and their corresponding labels
    directions = {(-1, 0): 'U', (1, 0): 'D', (0, -1): 'L', (0, 1): 'R'}
    
    found_path = None
    while queue:
        (r, c), path = queue.popleft()

        if (r, c) == goal_pos:
            found_path = path
            break

        # Explore neighbors
        for (dr, dc), move_char in directions.items():
            nr, nc = r + dr, c + dc
            
            # Check grid boundaries
            if not (0 <= nr < rows and 0 <= nc < cols):
                continue
            
            # Check if the tile has been visited or is a wall
            if (nr, nc) in visited:
                continue
            if grid[nr][nc] not in walkable_chars:
                continue

            # Add valid neighbor to the queue and mark as visited
            visited.add((nr, nc))
            new_path = path + [move_char]
            queue.append(((nr, nc), new_path))
    
    # Process and print the result
    if found_path:
        # Simplify the path by collapsing consecutive identical moves
        if not found_path:
            simplified_path_str = ""
        else:
            simplified_path = [found_path[0]]
            for move in found_path[1:]:
                if move != simplified_path[-1]:
                    simplified_path.append(move)
            simplified_path_str = "".join(simplified_path)
        
        print(f"The least dangerous path is: {simplified_path_str}")
        # The final answer in the required format
        print(f"<<<{simplified_path_str}>>>")
    else:
        print("No path to the gold could be found.")

solve_dungeon_path()