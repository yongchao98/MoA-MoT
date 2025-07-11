import collections

def solve_dungeon_path():
    """
    Solves the dungeon path problem by finding the shortest, safest route
    from the adventurer '@' to the gold 'g', avoiding the dragon 'D'.
    """

    # Raw text representation of the dungeon map
    raw_grid_text = r"""
    / &  & - & - & - & - & - & - &  &  &  &  &  &  &  &  &  &  &  &  &  &  &  &  &  &  &  &  &  &  \\
    / &  & | & . & . & . & . & | &  &  &  &  &  &  & \# & \# & \# & \# & \# & \# & \# & \# & \# & \# & \# & \# &  &  &  &  \\
    / &  & | & . & . & . & . & | &  &  &  &  &  &  & \# &  &  &  &  &  &  &  &  &  &  & \# &  &  &  &  \\
    / &  & | & . & g & . & . & + & \# & \# & \# & \# & \# & \# & \# & \# &  &  &  &  &  &  &  &  &  & @ &  &  &  &  \\
    / &  & | & . & . & . & . & | &  &  &  &  &  &  &  & \# &  &  &  &  &  &  & - & - & - & + & - & - & - &  \\
    / &  & - & - & - & - & - & - &  &  &  &  &  &  &  & \# &  &  &  &  &  &  & | & . & . & . & . & . & | \\
    / &  &  &  &  &  &  &  &  &  &  &  &  &  &  & \# &  &  &  &  &  &  & | & . & ! & . & . & . & | \\
    / &  &  &  &  &  &  &  &  &  &  &  &  &  &  & \# &  &  &  &  &  &  & | & . & . & . & . & . & | \\
    / &  &  &  &  &  &  &  &  &  &  &  &  &  &  & \# &  &  &  &  &  &  & | & . & . & . & . & . & | \\
    / &  &  &  & - & - & - & - &  &  &  &  &  &  &  & \# &  &  &  &  &  &  & | & . & . & . & . & . & | \\
    / &  &  &  & | & . & . & | &  &  &  &  &  &  &  & \# & \# & \# & \# & \# & \# & \# & + & . & . & D & . & . & | \\
    / &  &  &  & | & < & . & + & \# & \# & \# &  &  &  &  & \# &  &  &  &  &  &  & | & . & . & . & . & . & | \\
    / &  &  &  & - & - & - & - &  &  & \# &  &  &  &  & \# &  &  &  &  &  &  & | & . & ? & . & . & . & | \\
    / &  &  &  &  &  &  &  &  &  & \# & \# & \# & \# & \# & \# &  &  &  &  &  &  & - & - & - & - & - & - & - \\
    """

    # Parse the raw text into a 2D grid
    lines = raw_grid_text.strip().split('\\\\\n')
    grid = []
    for line in lines:
        line = line.strip()
        if not line:
            continue
        # Split line by '&' and replace empty parts with spaces
        cells = [cell.strip() if cell.strip() else ' ' for cell in line.split('&')]
        grid.append(cells)

    # Find the positions of the adventurer, gold, and dragon
    start_pos, goal_pos, dragon_pos = None, None, None
    for r, row in enumerate(grid):
        for c, char in enumerate(row):
            if char == '@':
                start_pos = (r, c)
            elif char == 'g':
                goal_pos = (r, c)
            elif char == 'D':
                dragon_pos = (r, c)

    if not start_pos or not goal_pos:
        print("Error: Could not find start ('@') or goal ('g').")
        return

    # Define obstacles and valid moves for BFS
    obstacles = {'#', '-', '|', '/'}
    moves = {'U': (-1, 0), 'D': (1, 0), 'L': (0, -1), 'R': (0, 1)}
    
    # Initialize the queue for BFS with (position, path_taken)
    queue = collections.deque([(start_pos, "")])
    visited = {start_pos}
    
    rows, cols = len(grid), len(grid[0])
    shortest_path = None

    # Perform Breadth-First Search
    while queue:
        (r, c), path = queue.popleft()

        # If the goal is reached, store the path and break
        if (r, c) == goal_pos:
            shortest_path = path
            break

        # Explore neighbors (Up, Down, Left, Right)
        for move_char, (dr, dc) in moves.items():
            nr, nc = r + dr, c + dc

            # Check if the new position is valid
            if 0 <= nr < rows and 0 <= nc < cols and \
               (nr, nc) not in visited and \
               grid[nr][nc] not in obstacles and \
               (nr, nc) != dragon_pos:
                
                visited.add((nr, nc))
                queue.append(((nr, nc), path + move_char))

    # Process and print the result
    if shortest_path is not None:
        # Compress the path sequence (e.g., "LLLL" becomes "L")
        if not shortest_path:
            # This case handles if start is the same as goal.
            print("Adventurer is already at the gold.")
        else:
            compressed_path = shortest_path[0]
            for i in range(1, len(shortest_path)):
                if shortest_path[i] != shortest_path[i-1]:
                    compressed_path += shortest_path[i]
            print(compressed_path)
    else:
        print("No safe path to the gold was found.")

solve_dungeon_path()