import heapq

def solve_sokoban_like_puzzle():
    """
    Solves the puzzle by finding the least dangerous path from '@' to 'g'.
    This function implements the A* search algorithm.
    """
    
    # 1. Parse the Map from the LaTeX bmatrix string
    latex_str = r"""
    / &  & - & - & - & - & - & - &  &  &  &  &  &  &  &  &  &  &  &  &  &  &  &  &  &  &  &  &  &  \\
    / &  & | & . & . & . & . & | &  &  &  &  &  &  & # & # & # & # & # & # & # & # & # & # & # & # &  &  &  &  \\
    / &  & | & . & . & . & . & | &  &  &  &  &  &  & # &  &  &  &  &  &  &  &  &  &  & # &  &  &  &  \\
    / &  & | & . & g & . & . & + & # & # & # & # & # & # & # & # &  &  &  &  &  &  &  &  &  & @ &  &  &  &  \\
    / &  & - & - & - & - & - & - &  &  &  &  &  &  &  & # &  &  &  &  &  &  & - & - & - & + & - & - & - &  \\
    / &  &  &  &  &  &  &  &  &  &  &  &  &  &  & # &  &  &  &  &  &  & | & . & . & . & . & . & | \\
    / &  &  &  &  &  &  &  &  &  &  &  &  &  &  & # &  &  &  &  &  &  & | & . & ! & . & . & . & | \\
    / &  &  &  &  &  &  &  &  &  &  &  &  &  &  & # &  &  &  &  &  &  & | & . & . & . & . & . & | \\
    / &  &  &  &  &  &  &  &  &  &  &  &  &  &  & # &  &  &  &  &  &  & | & . & . & . & . & . & | \\
    / &  &  &  & - & - & - & - &  &  &  &  &  &  &  & # &  &  &  &  &  &  & | & . & . & . & . & . & | \\
    / &  &  &  & | & . & . & | &  &  &  &  &  &  &  & # & # & # & # & # & # & # & + & . & . & D & . & . & | \\
    / &  &  &  & | & < & . & + & # & # & # &  &  &  &  & # &  &  &  &  &  &  & | & . & . & . & . & . & | \\
    / &  &  &  & - & - & - & - &  &  & # &  &  &  &  & # &  &  &  &  &  &  & | & . & ? & . & . & . & | \\
    / &  &  &  &  &  &  &  &  &  & # & # & # & # & # & # &  &  &  &  &  &  & - & - & - & - & - & - & - \\
    """
    
    lines = latex_str.strip().replace('\\\\', '\n').split('\n')
    grid = []
    for line in lines:
        if not line.strip(): continue
        cells = [c.strip() for c in line.split('&')]
        grid.append(cells)

    # 2. Define walkability and find key locations
    walkable_chars = {'.', '#', '+', 'g', '@', 'D', '!', '?', '<', ''}
    start_pos, goal_pos, dragon_pos = None, None, None
    for r, row in enumerate(grid):
        for c, char in enumerate(row):
            if char == '@':
                start_pos = (r, c)
            elif char == 'g':
                goal_pos = (r, c)
            elif char == 'D':
                dragon_pos = (r, c)

    rows, cols = len(grid), len(grid[0])

    def is_walkable(r, c):
        if not (0 <= r < rows and 0 <= c < cols):
            return False
        return grid[r][c] in walkable_chars

    # 3. Implement A*
    def manhattan_distance(p1, p2):
        return abs(p1[0] - p2[0]) + abs(p1[1] - p2[1])

    # A large factor to make dragon avoidance a high priority
    DANGER_FACTOR = 5000 
    
    # Priority queue stores: (priority, g_cost, path, current_pos)
    pq = [(0, 0, [start_pos], start_pos)]
    visited = {start_pos}

    final_path = None

    while pq:
        _, g_cost, path, pos = heapq.heappop(pq)
        
        if pos == goal_pos:
            final_path = path
            break

        r, c = pos
        moves = [(r-1, c), (r+1, c), (r, c-1), (r, c+1)] # U, D, L, R

        for next_pos in moves:
            if next_pos in visited or not is_walkable(next_pos[0], next_pos[1]):
                continue
            
            visited.add(next_pos)
            new_g_cost = g_cost + 1
            h_cost = manhattan_distance(next_pos, goal_pos)
            
            danger_cost = 0
            if dragon_pos:
                dist_to_dragon = manhattan_distance(next_pos, dragon_pos)
                danger_cost = DANGER_FACTOR / (dist_to_dragon + 1) # +1 to avoid division by zero
                
            priority = new_g_cost + h_cost + danger_cost
            new_path = path + [next_pos]
            heapq.heappush(pq, (priority, new_g_cost, new_path, next_pos))

    # 4. Convert path to simplified directions
    if not final_path:
        print("No path found.")
        return

    directions = []
    for i in range(1, len(final_path)):
        r1, c1 = final_path[i-1]
        r2, c2 = final_path[i]
        if r2 < r1: directions.append('U')
        elif r2 > r1: directions.append('D')
        elif c2 < c1: directions.append('L')
        elif c2 > c1: directions.append('R')
        
    if not directions:
        print("Adventurer is already at the gold.")
        return

    simplified_directions = [directions[0]]
    for i in range(1, len(directions)):
        if directions[i] != simplified_directions[-1]:
            simplified_directions.append(directions[i])
            
    result = "".join(simplified_directions)
    print(f"<<<{result}>>>")

solve_sokoban_like_puzzle()