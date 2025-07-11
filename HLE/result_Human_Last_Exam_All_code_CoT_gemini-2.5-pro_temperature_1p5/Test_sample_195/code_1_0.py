import heapq

def solve_sokoban_path():
    """
    Solves the maze problem by finding the least dangerous path.
    The "least dangerous" path is interpreted as the path that stays as far away
    from the dragon 'D' as possible. This is modeled as a shortest path problem
    on a weighted graph using Dijkstra's algorithm, where tiles near the
    dragon have a very high travel cost.
    """
    grid = [
        "/                                ",
        "/  | . . . . . |       # # # # # # # # # # # # # # # #",
        "/  | . . . . . |       #            #            ",
        "/  | . g . . . + # # # # # # # # # # # @            ",
        "/  | . . . . . |       #           - - - + - - - ",
        "/  - - - - - - -       #           | . . . . . | ",
        "/                        #           | . ! . . . | ",
        "/                        #           | . . . . . | ",
        "/                        #           | . . . . . | ",
        "/          - - - - - #           | . . . . . | ",
        "/          | . . | # # # # #     # + . . D . . | ",
        "/          | < . + # # # # # #   # | . . . . . | ",
        "/          - - - - -   #   #           | . ? . . . | ",
        "/                      # # # # # # # # # - - - - - - - - ",
    ]

    # Normalize grid dimensions
    max_len = max(len(row) for row in grid)
    grid = [row.ljust(max_len) for row in grid]

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

    if not all([start_pos, goal_pos, dragon_pos]):
        print("Error: Could not find start, goal, or dragon in the grid.")
        return

    walls = {'#', '|', '-'}
    
    # Priority queue for Dijkstra's: (cost, r, c)
    pq = [(0, start_pos)]
    
    # Dictionaries to store costs and path
    cost_so_far = {start_pos: 0}
    came_from = {start_pos: None}

    while pq:
        current_cost, current_pos = heapq.heappop(pq)

        if current_pos == goal_pos:
            break

        if current_cost > cost_so_far.get(current_pos, float('inf')):
            continue

        r, c = current_pos
        for dr, dc in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
            next_pos = (r + dr, c + dc)
            nr, nc = next_pos

            if not (0 <= nr < rows and 0 <= nc < cols):
                continue
            
            char_at_next = grid[nr][nc]
            if char_at_next in walls or char_at_next == 'D':
                continue

            # Calculate danger cost
            dist_to_dragon = max(abs(nr - dragon_pos[0]), abs(nc - dragon_pos[1]))
            
            move_cost = 1
            if dist_to_dragon <= 1: # Adjacent to dragon
                move_cost = 100
            elif dist_to_dragon <= 3: # Near dragon
                 move_cost = 10

            new_cost = current_cost + move_cost

            if new_cost < cost_so_far.get(next_pos, float('inf')):
                cost_so_far[next_pos] = new_cost
                came_from[next_pos] = current_pos
                heapq.heappush(pq, (new_cost, next_pos))

    # Reconstruct path
    path = []
    curr = goal_pos
    if curr not in came_from:
        print("No path found.")
        return

    while curr is not None:
        path.append(curr)
        curr = came_from[curr]
    path.reverse()

    # Convert path to directions
    directions = []
    for i in range(len(path) - 1):
        r1, c1 = path[i]
        r2, c2 = path[i+1]
        if r2 > r1: directions.append('D')
        elif r2 < r1: directions.append('U')
        elif c2 > c1: directions.append('R')
        elif c2 < c1: directions.append('L')

    # Simplify directions
    if not directions:
        print("No moves needed.")
        return
        
    simplified_directions = [directions[0]]
    for dir in directions[1:]:
        if dir != simplified_directions[-1]:
            simplified_directions.append(dir)

    print("".join(simplified_directions))

solve_sokoban_path()
<<<DLUL>>>