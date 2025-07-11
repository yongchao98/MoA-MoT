import heapq

def solve():
    """
    Solves the dungeon pathfinding problem.
    Finds the least dangerous path from the adventurer '@' to the gold '$g',
    avoiding the dragon 'D'.
    """
    map_text = r"""
    / &  & - & - & - & - & - & - &  &  &  &  &  &  &  &  &  &  &  &  &  &  &  &  &  &  &  &  &  &  \\
    / &  & | & . & . & . & . & | &  &  &  &  &  &  & # & # & # & # & # & # & # & # & # & # & # & # &  &  &  &  \\
    / &  & | & . & . & . & . & | &  &  &  &  &  &  & # &  &  &  &  &  &  &  &  &  &  & # &  &  &  &  \\
    / &  & | & . & g & . & . & + & # & # & # & # & # & # & # & # &  &  &  &  &  &  &  &  &  & @ &  &  &  &  \\
    / &  & | & . & . & . & . & | &  &  &  &  &  &  &  & # &  &  &  &  &  &  & - & - & - & + & - & - & - &  \\
    / &  & - & - & - & - & - & - &  &  &  &  &  &  &  & # &  &  &  &  &  &  & | & . & . & . & . & . & | \\
    / &  &  &  &  &  &  &  &  &  &  &  &  &  &  &  & # &  &  &  &  &  &  & | & . & ! & . & . & . & | \\
    / &  &  &  &  &  &  &  &  &  &  &  &  &  &  &  & # &  &  &  &  &  &  & | & . & . & . & . & . & | \\
    / &  &  &  &  &  &  &  &  &  &  &  &  &  &  &  & # &  &  &  &  &  &  & | & . & . & . & . & . & | \\
    / &  &  &  & - & - & - & - &  &  &  &  &  &  &  &  & # &  &  &  &  &  &  & | & . & . & . & . & . & | \\
    / &  &  &  & | & . & . & | &  &  &  &  &  &  &  &  & # & # & # & # & # & # & + & . & . & D & . & . & | \\
    / &  &  &  & | & < & . & + & # & # & # &  &  &  &  &  & # &  &  &  &  &  &  & | & . & . & . & . & . & | \\
    / &  &  &  & - & - & - & - &  &  & # &  &  &  &  &  & # &  &  &  &  &  &  & | & . & ? & . & . & . & | \\
    / &  &  &  &  &  &  &  &  &  & # & # & # & # & # & # &  &  &  &  &  &  &  & - & - & - & - & - & - & - \\
    """

    def parse_map(text):
        """Parses the ASCII map into a 2D grid."""
        grid = []
        lines = text.strip().split('\\')
        for line in lines:
            cells = [c.strip() for c in line.strip().replace('/', '').split('&')]
            grid.append(cells)
        # Normalize grid by padding rows to the same length
        max_len = 0
        if grid:
            max_len = max(len(r) for r in grid)
        for row in grid:
            row.extend([' '] * (max_len - len(row)))
        return grid

    grid = parse_map(map_text)
    
    start_pos, goal_pos, dragon_pos = None, None, None
    for r, row in enumerate(grid):
        for c, char in enumerate(row):
            if char == '@':
                start_pos = (r, c)
            elif char == 'g':
                goal_pos = (r, c)
            elif char == 'D':
                dragon_pos = (r, c)

    def heuristic(a, b):
        """Calculates Manhattan distance heuristic."""
        return abs(a[0] - b[0]) + abs(a[1] - b[1])

    def get_cost(pos, dragon_pos):
        """Calculates movement cost for a tile, penalizing proximity to the dragon."""
        char = grid[pos[0]][pos[1]]
        if char == '#' or char == 'D':
            return float('inf')
        
        danger_penalty = 0
        if dragon_pos:
            dist = max(abs(pos[0] - dragon_pos[0]), abs(pos[1] - dragon_pos[1]))
            if dist <= 2:
                danger_penalty = 100 # High cost for being near the dragon
        return 1 + danger_penalty

    def a_star_search(start, goal):
        """Performs A* search to find the least dangerous path."""
        open_set = [(0, start)]  # (f_cost, pos)
        came_from = {}
        g_score = {start: 0}
        f_score = {start: heuristic(start, goal)}

        while open_set:
            _, current = heapq.heappop(open_set)

            if current == goal:
                return reconstruct_path(came_from, current)

            for dr, dc in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
                neighbor = (current[0] + dr, current[1] + dc)
                
                if not (0 <= neighbor[0] < len(grid) and 0 <= neighbor[1] < len(grid[0])):
                    continue

                cost = get_cost(neighbor, dragon_pos)
                if cost == float('inf'):
                    continue

                tentative_g_score = g_score[current] + cost
                
                if neighbor not in g_score or tentative_g_score < g_score[neighbor]:
                    came_from[neighbor] = current
                    g_score[neighbor] = tentative_g_score
                    f_score[neighbor] = tentative_g_score + heuristic(neighbor, goal)
                    heapq.heappush(open_set, (f_score[neighbor], neighbor))
        return None

    def reconstruct_path(came_from, current):
        """Reconstructs the path from the came_from map."""
        path = [current]
        while current in came_from:
            current = came_from[current]
            path.append(current)
        return path[::-1]

    def simplify_path(path):
        """Converts a path of coordinates into a simplified sequence of moves."""
        if not path or len(path) < 2:
            return ""

        moves = []
        for i in range(len(path) - 1):
            r1, c1 = path[i]
            r2, c2 = path[i+1]
            if r2 > r1: moves.append('D')
            elif r2 < r1: moves.append('U')
            elif c2 > c1: moves.append('R')
            elif c2 < c1: moves.append('L')

        if not moves:
            return ""

        simplified_moves = [moves[0]]
        for i in range(1, len(moves)):
            if moves[i] != moves[i-1]:
                simplified_moves.append(moves[i])
        
        return "".join(simplified_moves)

    path = a_star_search(start_pos, goal_pos)
    
    if path:
        simplified_path = simplify_path(path)
        print(simplified_path)
    else:
        print("No path found.")

solve()
<<<DLDLUL>>>