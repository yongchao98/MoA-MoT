import heapq
import math

def solve_maze():
    """
    Analyzes the map to find the least dangerous path from '@' to 'g'.
    The pathfinding considers distance to the dragon 'D' as the primary danger factor.
    """
    # 1. Map Representation: The map is transcribed into a 2D grid.
    grid_str = [
        "                               ",
        "   |.......|  ###########      ",
        "   |.......|  #         #      ",
        "   |..g..+###########  @        ",
        "   |.......|  #   ---+-        ",
        "   ---------  #   |.....|      ",
        "              #   |.!...|      ",
        "              #   |.....|      ",
        "              #   |.....|      ",
        "        ----  #   |.....|      ",
        "        |..|  #######+..D..|   ",
        "        |<.+###  #   |.....|   ",
        "        ----  #  #   |.?...|   ",
        "              ###########      "
    ]
    # Ensure the grid is rectangular by padding with spaces.
    max_len = max(len(row) for row in grid_str)
    grid = [row.ljust(max_len) for row in grid_str]
    height, width = len(grid), len(grid[0])

    # 2. Identify Key Points
    start_pos, goal_pos, dragon_pos = None, None, None
    for r, row in enumerate(grid):
        for c, char in enumerate(row):
            if char == '@': start_pos = (r, c)
            elif char == 'g': goal_pos = (r, c)
            elif char == 'D': dragon_pos = (r, c)

    if not all((start_pos, goal_pos, dragon_pos)):
        print("Map error: Could not locate '@', 'g', or 'D'.")
        return

    # 3. Pathfinding Algorithm (Dijkstra's)
    # The priority queue stores tuples of (cost, row, col, path_history).
    pq = [(0, start_pos[0], start_pos[1], [])]
    # visited keeps track of the minimum cost found so far to reach a cell.
    visited = {}

    walls = {'-', '|', ' '}

    while pq:
        cost, r, c, path = heapq.heappop(pq)

        if (r, c) in visited and visited[(r, c)] <= cost:
            continue
        visited[(r, c)] = cost

        new_path = path + [(r, c)]

        if (r, c) == goal_pos:
            # 5. Path Simplification and 6. Final Output
            convert_path_to_directions(new_path)
            return

        # Explore neighbors (Up, Down, Left, Right)
        for dr, dc in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
            nr, nc = r + dr, c + dc

            if not (0 <= nr < height and 0 <= nc < width):
                continue

            tile = grid[nr][nc]
            if tile in walls:
                continue

            # 4. Cost Calculation
            # Base cost: Unlit '#' is slightly more expensive.
            move_cost = 2 if tile == '#' else 1

            # Danger cost: A large penalty for proximity to the dragon.
            dist_sq = (nr - dragon_pos[0])**2 + (nc - dragon_pos[1])**2
            danger_penalty = 1000 / (dist_sq + 1e-9)  # Add epsilon to avoid division by zero.

            new_cost = cost + move_cost + danger_penalty

            if (nr, nc) not in visited or new_cost < visited[(nr, nc)]:
                heapq.heappush(pq, (new_cost, nr, nc, new_path))

    print("No path to the gold could be found.")

def convert_path_to_directions(path_coords):
    """
    Converts a list of coordinates into a simplified string of moves and prints it.
    """
    if len(path_coords) < 2:
        return

    moves = []
    for i in range(len(path_coords) - 1):
        r1, c1 = path_coords[i]
        r2, c2 = path_coords[i+1]
        if r2 > r1: moves.append('D')
        elif r2 < r1: moves.append('U')
        elif c2 > c1: moves.append('R')
        elif c2 < c1: moves.append('L')

    if not moves:
        return

    # Condense sequential moves like L, L, L into a single L.
    simplified_moves = [moves[0]]
    for move in moves[1:]:
        if move != simplified_moves[-1]:
            simplified_moves.append(move)

    print("".join(simplified_moves))

# Execute the solver.
solve_maze()