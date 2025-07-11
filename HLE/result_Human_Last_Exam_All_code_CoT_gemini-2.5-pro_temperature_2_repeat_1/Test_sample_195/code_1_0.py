import collections

def solve_dungeon_path():
    """
    Solves the dungeon pathfinding problem by finding the least dangerous
    path from the adventurer '@' to the gold 'g'.
    """
    raw_map_str = r"""
/ &  & - & - & - & - & - & - &  &  &  &  &  &  &  &  &  &  &  &  &  &  &  &  &  &  &  &  &  &  \\
/ &  & | & . & . & . & . & | &  &  &  &  &  &  & # & # & # & # & # & # & # & # & # & # & # & # &  &  &  &  \\
/ &  & | & . & . & . & . & | &  &  &  &  &  &  & # &  &  &  &  &  &  &  &  &  &  & # &  &  &  &  \\
/ &  & | & . & g & . & . & + & # & # & # & # & # & # & # & # &  &  &  &  &  &  &  &  &  & @ &  &  &  &  \\
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
    for line in lines:
        line = line.strip()
        if not line:
            continue
        # Each cell is separated by '&'. We skip the first column ('/').
        parts = [p.strip() for p in line.split('&')]
        grid.append(parts[1:])

    rows, cols = len(grid), len(grid[0])
    start_pos, goal_pos, dragon_pos = None, None, None

    for r in range(rows):
        for c in range(cols):
            if grid[r][c] == '@':
                start_pos = (r, c)
            elif grid[r][c] == 'g':
                goal_pos = (r, c)
            elif grid[r][c] == 'D':
                dragon_pos = (r, c)
    
    walkable_chars = ['.', '#', '+', 'g', '!', '?', '<', '@']
    
    # Queue for BFS: stores tuples of (position, path_string)
    queue = collections.deque([(start_pos, "")])
    visited = {start_pos}

    path_to_goal = None

    while queue:
        (r, c), path = queue.popleft()

        if (r, c) == goal_pos:
            path_to_goal = path
            break

        # Explore neighbors: Right, Left, Down, Up
        for dr, dc, move in [(0, 1, 'R'), (0, -1, 'L'), (1, 0, 'D'), (-1, 0, 'U')]:
            nr, nc = r + dr, c + dc

            if 0 <= nr < rows and 0 <= nc < cols and (nr, nc) not in visited:
                char_at_new_pos = grid[nr][nc]
                
                # The "least dangerous" path avoids the dragon
                if (nr, nc) == dragon_pos:
                    continue
                
                # Check if the tile is walkable
                if char_at_new_pos in walkable_chars:
                    visited.add((nr, nc))
                    queue.append(((nr, nc), path + move))

    if path_to_goal:
        # Condense the path, e.g., 'DDD' becomes 'D'
        if not path_to_goal:
            print("Path is empty.")
            return

        condensed_path = path_to_goal[0]
        for i in range(1, len(path_to_goal)):
            if path_to_goal[i] != path_to_goal[i-1]:
                condensed_path += path_to_goal[i]
        
        print(condensed_path)
    else:
        print("No path found.")

solve_dungeon_path()