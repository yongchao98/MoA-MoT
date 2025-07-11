import heapq

def find_least_dangerous_path():
    """
    This function solves the maze problem by:
    1.  Parsing the map into a 2D grid.
    2.  Identifying the start ('@'), goal ('g'), and danger ('D').
    3.  Defining a "danger cost" for moving onto tiles. The primary danger is the dragon's room.
    4.  Using Dijkstra's algorithm to find the shortest path from start to goal that avoids high-danger areas.
    5.  Converting the path into a compressed sequence of moves (U, D, L, R).
    """

    # Step 1: Parse the map into a 2D grid.
    # The map is given in a format similar to a LaTeX table, with '&' as a column separator.
    # We will split each line by '&' to preserve the column structure.
    raw_map_data = [
        r'/ &  & - & - & - & - & - & - &  &  &  &  &  &  &  &  &  &  &  &  &  &  &  &  &  &  &  &  &  &  \\',
        r'/ &  & | & . & . & . & . & | &  &  &  &  &  &  & # & # & # & # & # & # & # & # & # & # & # & # &  &  &  &  \\',
        r'/ &  & | & . & . & . & . & | &  &  &  &  &  &  & # &  &  &  &  &  &  &  &  &  &  & # &  &  &  &  \\',
        r'/ &  & | & . & g & . & . & + & # & # & # & # & # & # & # & # &  &  &  &  &  &  &  &  &  & @ &  &  &  &  \\',
        r'/ &  & | & . & . & . & . & | &  &  &  &  &  &  &  & # &  &  &  &  &  &  & - & - & - & + & - & - & - &  \\',
        r'/ &  & - & - & - & - & - & - &  &  &  &  &  &  &  & # &  &  &  &  &  &  & | & . & . & . & . & . & | \\',
        r'/ &  &  &  &  &  &  &  &  &  &  &  &  &  &  & # &  &  &  &  &  &  & | & . & ! & . & . & . & | \\',
        r'/ &  &  &  &  &  &  &  &  &  &  &  &  &  &  & # &  &  &  &  &  &  & | & . & . & . & . & . & | \\',
        r'/ &  &  &  &  &  &  &  &  &  &  &  &  &  &  & # &  &  &  &  &  &  & | & . & . & . & . & . & | \\',
        r'/ &  &  &  & - & - & - & - &  &  &  &  &  &  &  & # &  &  &  &  &  &  & | & . & . & . & . & . & | \\',
        r'/ &  &  &  & | & . & . & | &  &  &  &  &  &  &  & # & # & # & # & # & # & # & + & . & . & D & . & . & | \\',
        r'/ &  &  &  & | & < & . & + & # & # & # &  &  &  &  & # &  &  &  &  &  &  & | & . & . & . & . & . & | \\',
        r'/ &  &  &  & - & - & - & - &  &  & # &  &  &  &  & # &  &  &  &  &  &  & | & . & ? & . & . & . & | \\',
        r'/ &  &  &  &  &  &  &  &  &  & # & # & # & # & # & # &  &  &  &  &  &  & - & - & - & - & - & - & - \\'
    ]
    grid = [[c.strip() or ' ' for c in line.strip().strip('\\').split('&')] for line in raw_map_data]

    # Find key locations
    rows, cols = len(grid), len(grid[0])
    start_pos, goal_pos = None, None
    for r in range(rows):
        for c in range(cols):
            if grid[r][c] == '@':
                start_pos = (r, c)
            elif grid[r][c] == 'g':
                goal_pos = (r, c)

    # Step 2: Define the cost of moving to a tile (our definition of "danger")
    def get_move_cost(r, c):
        char = grid[r][c]
        # Walls are impassable
        if char in ['|', '-', '/']:
            return float('inf')
        # The dragon's room is extremely dangerous. We define its area heuristically.
        # This area is below the main level and to the right of the central hallway.
        if r > 4 and c > 22:
            return 1000000
        # All other walkable tiles have a base cost of 1.
        return 1

    # Step 3: Use Dijkstra's algorithm to find the lowest-cost path
    # Priority queue stores: (total_cost, current_position)
    pq = [(0, start_pos)]
    # 'came_from' will store the path for later reconstruction
    came_from = {start_pos: None}
    # 'cost_so_far' stores the minimum cost found so far to reach a node
    cost_so_far = {start_pos: 0}

    path_found = False
    while pq:
        cost, current_pos = heapq.heappop(pq)

        if current_pos == goal_pos:
            path_found = True
            break
        
        # Check all four neighbors (Up, Down, Left, Right)
        r, c = current_pos
        for dr, dc in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
            nr, nc = r + dr, c + dc

            # Check if neighbor is within grid bounds
            if not (0 <= nr < rows and 0 <= nc < cols):
                continue
            
            # Calculate cost to move to neighbor
            move_cost = get_move_cost(nr, nc)
            new_cost = cost_so_far[current_pos] + move_cost
            
            # If the neighbor is not a wall and we found a better path
            if (nr, nc) not in cost_so_far or new_cost < cost_so_far[(nr, nc)]:
                cost_so_far[(nr, nc)] = new_cost
                priority = new_cost
                heapq.heappush(pq, (priority, (nr, nc)))
                came_from[(nr, nc)] = current_pos

    if not path_found:
        print("No path to the gold could be found.")
        return

    # Step 4: Reconstruct the path from the 'came_from' map
    path = []
    curr = goal_pos
    while curr is not None:
        path.append(curr)
        curr = came_from.get(curr)
    path.reverse()

    # Step 5: Convert the coordinate path to a sequence of U, D, L, R moves
    moves = []
    for i in range(len(path) - 1):
        r1, c1 = path[i]
        r2, c2 = path[i+1]
        if r2 > r1: moves.append('D')
        elif r2 < r1: moves.append('U')
        elif c2 > c1: moves.append('R')
        elif c2 < c1: moves.append('L')

    # Compress the move sequence (e.g., LLLLUU -> LU)
    if not moves:
        print("Already at the destination.")
        return
        
    compressed_moves = [moves[0]]
    for move in moves[1:]:
        if move != compressed_moves[-1]:
            compressed_moves.append(move)
    
    final_answer = "".join(compressed_moves)
    print(final_answer)

find_least_dangerous_path()