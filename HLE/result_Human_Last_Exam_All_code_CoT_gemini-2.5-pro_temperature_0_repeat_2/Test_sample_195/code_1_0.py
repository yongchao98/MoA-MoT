import heapq

def solve_maze():
    """
    This function solves the maze by finding the least dangerous path from the adventurer '@'
    to the gold 'g', avoiding the dragon 'D'. It uses Dijkstra's algorithm on a grid
    representation of the map, where the cost of traversing a tile is increased
    by its proximity to the dragon.
    """
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

    # 1. Parse the map into a 2D grid
    lines = grid_str.strip().split('\\')
    grid = []
    max_cols = 0
    parsed_lines = []
    for line in lines:
        line = line.strip()
        if not line: continue
        cells = [c.strip() for c in line.split('&')]
        parsed_lines.append(cells)
        max_cols = max(max_cols, len(cells))

    for row in parsed_lines:
        row.extend([' '] * (max_cols - len(row)))
        grid.append(row)

    rows, cols = len(grid), len(grid[0])
    start_pos, goal_pos, dragon_pos = None, None, None
    for r in range(rows):
        for c in range(cols):
            if grid[r][c] == '@': start_pos = (r, c)
            elif grid[r][c] == 'g': goal_pos = (r, c)
            elif grid[r][c] == 'D': dragon_pos = (r, c)

    # 2. Define pathfinding parameters
    def is_obstacle(r, c):
        if not (0 <= r < rows and 0 <= c < cols):
            return True
        return grid[r][c] in ['#', '/']

    def get_cost(pos):
        dist_to_dragon = abs(pos[0] - dragon_pos[0]) + abs(pos[1] - dragon_pos[1])
        danger_cost = 0
        if dist_to_dragon > 0:
            danger_cost = 1000 / (dist_to_dragon ** 2)
        else:
            danger_cost = float('inf')
        return 1 + danger_cost

    # 3. Run Dijkstra's algorithm
    pq = [(0, start_pos, [])]
    visited_costs = {start_pos: 0}
    path_found = None

    while pq:
        cost, current_pos, path = heapq.heappop(pq)

        if cost > visited_costs.get(current_pos, float('inf')):
            continue

        new_path = path + [current_pos]

        if current_pos == goal_pos:
            path_found = new_path
            break

        r, c = current_pos
        for dr, dc in [(-1, 0), (1, 0), (0, -1), (0, 1)]:  # U, D, L, R
            next_pos = (r + dr, c + dc)
            
            if not is_obstacle(next_pos[0], next_pos[1]):
                move_cost = get_cost(next_pos)
                new_cost = cost + move_cost
                
                if new_cost < visited_costs.get(next_pos, float('inf')):
                    visited_costs[next_pos] = new_cost
                    heapq.heappush(pq, (new_cost, next_pos, new_path))

    # 4. Convert path to simplified directions
    if not path_found:
        print("No path found.")
        return

    directions = []
    for i in range(len(path_found) - 1):
        r1, c1 = path_found[i]
        r2, c2 = path_found[i+1]
        if r2 > r1: directions.append('D')
        elif r2 < r1: directions.append('U')
        elif c2 > c1: directions.append('R')
        elif c2 < c1: directions.append('L')

    if not directions:
        print("")
        return
        
    simplified_path = [directions[0]]
    for i in range(1, len(directions)):
        if directions[i] != simplified_path[-1]:
            simplified_path.append(directions[i])
    
    final_path_str = "".join(simplified_path)
    print(final_path_str)
    print(f"<<<{final_path_str}>>>")

solve_maze()