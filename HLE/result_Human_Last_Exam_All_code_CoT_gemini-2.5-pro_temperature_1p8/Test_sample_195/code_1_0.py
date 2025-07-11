import collections

def solve_maze():
    """
    Solves the maze to find the least dangerous path to the gold.
    """
    # Step 1: Parse the map from its string representation
    raw_map_str = r"""
/ &  & - & - & - & - & - & - &  &  &  &  &  &  &  &  &  &  &  &  &  &  &  &  &  &  &  &  &  &  \\
/ &  & | & . & . & . & . & | &  &  &  &  &  &  & # & # & # & # & # & # & # & # & # & # & # & # &  &  &  &  \\
/ &  & | & . & . & . & . & | &  &  &  &  &  &  & # &  &  &  &  &  &  &  &  &  &  & # &  &  &  &  \\
/ &  & | & . & g & . & . & + & # & # & # & # & # & # & # &  &  &  &  &  &  &  &  &  & @ &  &  &  &  \\
/ &  & | & . & . & . & . & | &  &  &  &  &  &  &  & # &  &  &  &  &  &  & - & - & - & + & - & - & - &  \\
/ &  & - & - & - & - & - & - &  &  &  &  &  &  &  & # &  &  &  &  &  &  & | & . & . & . & . & . & | \\
/ &  &  &  &  &  &  &  &  &  &  &  &  &  &  & # &  &  &  &  &  &  & | & . & ! & . & . & . & | \\
/ &  &  &  &  &  &  &  &  &  &  &  &  &  &  & # &  &  &  &  &  &  & | & . & . & . & . & . & | \\
/ &  &  &  &  &  &  &  &  &  &  &  &  &  &  & # &  &  &  &  &  &  & | & . & . & . & . & . & | \\
/ &  &  &  & - & - & - & - &  &  &  &  &  &  &  & # &  &  &  &  &  &  & | & . & . & . & . & . & | \\
/ &  &  &  & | & . & . & | &  &  &  &  &  &  &  & # & # & # & # & # & # & # & + & . & . & D & . & . & | \\
/ &  &  &  & | & < & . & + & # & # & # &  &  &  &  & # &  &  &  &  &  &  & | & . & . & . & . & . & | \\
/ &  &  &  & - & - & - & - &  &  & # &  &  &  &  & # &  &  &  &  &  &  & | & . & ? & . & . & . & | \\
/ &  &  &  &  &  &  &  &  &  & # & # & # & # & # & # &  &  &  &  &  &  & - & - & - & - & - & - & - \\
    """
    lines = raw_map_str.strip().split('\\\\')
    grid = []
    max_width = 0
    for line in lines:
        # By replacing '&' with nothing, we preserve the spacing and alignment
        clean_line = line.strip().lstrip('/').replace('&', '')
        grid.append(list(clean_line))
        if len(clean_line) > max_width:
            max_width = len(clean_line)

    # Normalize grid to be rectangular by padding with spaces
    for i in range(len(grid)):
        row_len = len(grid[i])
        if row_len < max_width:
            grid[i].extend([' '] * (max_width - row_len))

    # Step 2: Find start, goal, and obstacles
    start_pos, goal_pos = None, None
    walls = {'-', '|', 'D'}
    height, width = len(grid), len(grid[0])
    for r in range(height):
        for c in range(width):
            if grid[r][c] == '@':
                start_pos = (r, c)
            elif grid[r][c] == 'g':
                goal_pos = (r, c)

    if not start_pos or not goal_pos:
        print("Start or goal not found.")
        return

    # Step 3: Find the shortest path using Breadth-First Search (BFS)
    queue = collections.deque([(start_pos, [])]) # (position, path_list)
    visited = {start_pos}
    path_found = None

    moves = {'U': (-1, 0), 'D': (1, 0), 'L': (0, -1), 'R': (0, 1)}

    while queue:
        (r, c), path = queue.popleft()

        if (r, c) == goal_pos:
            path_found = path
            break

        for move_char, (dr, dc) in moves.items():
            nr, nc = r + dr, c + dc

            if 0 <= nr < height and 0 <= nc < width and (nr, nc) not in visited:
                if grid[nr][nc] not in walls:
                    visited.add((nr, nc))
                    new_path = path + [move_char]
                    queue.append(((nr, nc), new_path))
    
    if not path_found:
        print("No path found to the gold.")
        return

    # Step 4: Simplify the path directions
    if not path_found:
        print("")
        return
        
    simplified_path = []
    if path_found:
        simplified_path.append(path_found[0])
        for move in path_found[1:]:
            if move != simplified_path[-1]:
                simplified_path.append(move)

    print("".join(simplified_path))

solve_maze()