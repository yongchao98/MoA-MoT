import collections

def solve_maze():
    """
    Solves the maze by finding the least dangerous path from '@' to 'g'.
    This is achieved by finding the shortest path using BFS, as the shortest
    path is also the safest one in this map layout.
    """
    # The map is manually transcribed into a grid.
    # '.' is used to represent any walkable space (lit areas, hallways).
    grid = [
        "//////////////////////////////////////",
        "/       -------         ############",
        "/       |.....|         #..........#",
        "/       |..g..+...................@#",
        "/       -------         #....---+-#",
        "/                       #....|.....|",
        "/                       #....|..!..|",
        "/                       #....|.....|",
        "/                       #....|.....|",
        "/         -----         #....|.....|",
        "/         |...|         #+#..D...|",
        "/         |<..+###       #....|.....|",
        "/         ----- ###       #....|.?...|",
        "/               #####       #----#####",
        "//////////////////////////////////////",
    ]

    # Find start (@) and goal (g) coordinates
    start_pos, goal_pos = None, None
    for r, row in enumerate(grid):
        for c, char in enumerate(row):
            if char == '@':
                start_pos = (r, c)
            elif char == 'g':
                goal_pos = (r, c)

    if not start_pos or not goal_pos:
        print("Start or goal not found in the map.")
        return

    # BFS implementation to find the shortest path
    queue = collections.deque([(start_pos, [])])
    visited = {start_pos}
    
    # U, D, L, R moves
    moves = {'U': (-1, 0), 'D': (1, 0), 'L': (0, -1), 'R': (0, 1)}
    
    # Define obstacles
    obstacles = {'#', '-', '|', '/'}
    
    path = None
    while queue:
        (r, c), current_path = queue.popleft()

        if (r, c) == goal_pos:
            path = current_path
            break

        for move_char, (dr, dc) in moves.items():
            nr, nc = r + dr, c + dc

            # Check if the new position is valid
            if 0 <= nr < len(grid) and 0 <= nc < len(grid[0]) and \
               (nr, nc) not in visited and \
               grid[nr][nc] not in obstacles:
                
                visited.add((nr, nc))
                new_path = current_path + [move_char]
                queue.append(((nr, nc), new_path))
    
    if not path:
        print("No path found.")
        return

    # Condense the path as per instructions
    # e.g., ['L', 'L', 'D'] becomes ['L', 'D']
    if not path:
        condensed_path_str = ""
    else:
        condensed_path = [path[0]]
        for move in path[1:]:
            if move != condensed_path[-1]:
                condensed_path.append(move)
        condensed_path_str = "".join(condensed_path)
    
    print(condensed_path_str)

solve_maze()
<<<L>>>