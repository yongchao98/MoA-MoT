import heapq
import math

def solve_dungeon():
    """
    Solves the dungeon puzzle by finding the least dangerous path from '@' to 'g'.
    """
    
    # 1. Parse the map into a grid
    grid_str = r"""
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
    
    grid = [list(line[2:-3]) for line in grid_str.strip().split('\n')]
    
    # Identify key locations
    start_pos, goal_pos, dragon_pos = None, None, None
    for r, row in enumerate(grid):
        for c, char in enumerate(row):
            if char == '@':
                start_pos = (r, c)
            elif char == 'g':
                goal_pos = (r, c)
            elif char == 'D':
                dragon_pos = (r, c)

    if not start_pos or not goal_pos or not dragon_pos:
        print("Error: Could not find start, goal, or dragon on the map.")
        return
        
    walls = {'#', '|', '-'}
    
    def get_cost(pos):
        """Calculate the danger cost of a tile based on proximity to the dragon."""
        dist_sq = (pos[0] - dragon_pos[0])**2 + (pos[1] - dragon_pos[1])**2
        # Danger cost is high near the dragon, falls off with square of distance
        # Add 1 to cost for movement, and a high penalty for being near the dragon
        return 1 + 500 / (1 + dist_sq)

    # 2. Use A* pathfinding algorithm
    # Heuristic: Manhattan distance to the goal
    def heuristic(a, b):
        return abs(a[0] - b[0]) + abs(a[1] - b[1])

    open_set = [(0, start_pos)]
    came_from = {}
    g_score = {start_pos: 0}
    
    while open_set:
        _, current = heapq.heappop(open_set)

        if current == goal_pos:
            break

        for dr, dc in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
            neighbor = (current[0] + dr, current[1] + dc)

            if not (0 <= neighbor[0] < len(grid) and 0 <= neighbor[1] < len(grid[0])):
                continue

            if grid[neighbor[0]][neighbor[1]] in walls:
                continue

            # Tentative g_score is the cost from start to current + danger cost of neighbor
            tentative_g_score = g_score[current] + get_cost(neighbor)

            if neighbor not in g_score or tentative_g_score < g_score.get(neighbor, float('inf')):
                came_from[neighbor] = current
                g_score[neighbor] = tentative_g_score
                f_score = tentative_g_score + heuristic(neighbor, goal_pos)
                heapq.heappush(open_set, (f_score, neighbor))
    
    # 3. Reconstruct and summarize the path
    path = []
    current = goal_pos
    while current in came_from:
        path.append(current)
        current = came_from[current]
    path.append(start_pos)
    path.reverse()
    
    if not path or len(path) < 2:
        print("No path found.")
        return

    # Identify major turns in the path
    turns = [path[0]]
    for i in range(1, len(path) - 1):
        r_prev, c_prev = path[i - 1]
        r_curr, c_curr = path[i]
        r_next, c_next = path[i + 1]

        dir1 = (r_curr - r_prev, c_curr - c_prev)
        dir2 = (r_next - r_curr, c_next - c_curr)

        if dir1 != dir2:
            turns.append(path[i])
    turns.append(path[-1])

    # Convert path segments into directions
    directions = []
    for i in range(len(turns) - 1):
        start_segment = turns[i]
        end_segment = turns[i+1]

        delta_r = end_segment[0] - start_segment[0]
        delta_c = end_segment[1] - start_segment[1]
        
        if abs(delta_r) > abs(delta_c): # Vertical movement
            directions.append('D' if delta_r > 0 else 'U')
        else: # Horizontal movement
            directions.append('R' if delta_c > 0 else 'L')

    print("".join(directions))

solve_dungeon()
<<<DLUL>>>