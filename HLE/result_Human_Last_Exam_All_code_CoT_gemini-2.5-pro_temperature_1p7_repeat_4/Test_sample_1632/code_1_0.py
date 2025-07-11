def solve_saw_count():
    """
    Calculates a(10), the number of 10-step self-avoiding walks
    on a Manhattan (square) lattice using a recursive approach.
    """

    # The target number of steps for the walk.
    n = 10

    # This recursive function counts all valid walks starting from the given path.
    # 'path' is the list of (x, y) coordinates making up the current walk.
    # 'visited' is a set of the same coordinates for fast O(1) lookup to check for self-intersections.
    def count_walks_recursive(path, visited):
        # The base case: if the path has n steps (i.e., n+1 points), we have found one valid walk.
        if len(path) - 1 == n:
            return 1

        count = 0
        # Get the last point in the path to start the next step from.
        x, y = path[-1]

        # Define the four possible moves: up, down, right, left.
        moves = [(0, 1), (0, -1), (1, 0), (-1, 0)]

        for dx, dy in moves:
            next_point = (x + dx, y + dy)

            # Check the "self-avoiding" constraint.
            if next_point not in visited:
                # If the point is new, take the step.
                path.append(next_point)
                visited.add(next_point)

                # Recursively count the walks from this new state.
                count += count_walks_recursive(path, visited)

                # Backtrack: undo the move so we can explore other directions from the previous point.
                path.pop()
                visited.remove(next_point)

        return count

    # All walks start at the origin (0, 0).
    initial_path = [(0, 0)]
    initial_visited = {(0, 0)}

    # Start the recursive counting process.
    result = count_walks_recursive(initial_path, initial_visited)

    # Per the instruction "output each number in the final equation",
    # we print the equation `a(n) = result`.
    print(f"a({n}) = {result}")

solve_saw_count()