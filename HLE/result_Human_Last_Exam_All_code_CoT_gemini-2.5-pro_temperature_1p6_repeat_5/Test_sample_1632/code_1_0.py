def solve_saw_count():
    """
    Calculates a(10), the number of 10-step self-avoiding walks on a Manhattan lattice.
    
    This function uses a recursive backtracking algorithm optimized by exploiting the
    symmetries of the square lattice to reduce the computation time.
    """

    def count_walks(steps_left, x, y, visited):
        """
        Recursively counts the number of self-avoiding walks from a given state.

        Args:
            steps_left: Number of steps remaining in the walk.
            x: The current x-coordinate.
            y: The current y-coordinate.
            visited: A set of (x, y) tuples representing visited points on the path.

        Returns:
            The total number of valid self-avoiding walks from the current state.
        """
        # Base case: If no steps are left, we have found one complete walk.
        if steps_left == 0:
            return 1

        total_walks = 0
        # Explore the four possible moves on a Manhattan lattice.
        for dx, dy in [(1, 0), (-1, 0), (0, 1), (0, -1)]:
            next_x, next_y = x + dx, y + dy
            
            # If the next point has not been visited, proceed with the walk.
            if (next_x, next_y) not in visited:
                # Add the new point to the visited set for the recursive call.
                visited.add((next_x, next_y))
                total_walks += count_walks(steps_left - 1, next_x, next_y, visited)
                # Backtrack: remove the point to explore other paths from the current state.
                visited.remove((next_x, next_y))

        return total_walks

    N = 10

    # The problem can be broken down using symmetry.
    # A walk's first step can be in 4 directions. Let's fix it to (0,0)->(1,0).
    # The second step can be "straight" to (2,0) or a "turn" to (1,1) or (1,-1).
    # The number of paths for a right turn and left turn are equal by symmetry.
    # So, a(N) = 4 * (C_straight + 2 * C_turn)
    # We need to count the remaining N-2 = 8 steps for each case.

    # Case 1: "Straight" start: (0,0) -> (1,0) -> (2,0)
    # We count 8-step walks starting from (2,0), avoiding the existing path.
    visited_straight = {(0, 0), (1, 0), (2, 0)}
    c_straight = count_walks(N - 2, 2, 0, visited_straight)

    # Case 2: "Turn" start: (0,0) -> (1,0) -> (1,1)
    # We count 8-step walks starting from (1,1), avoiding the existing path.
    visited_turn = {(0, 0), (1, 0), (1, 1)}
    c_turn = count_walks(N - 2, 1, 1, visited_turn)
    
    # Calculate the total using the symmetry formula.
    total = 4 * (c_straight + 2 * c_turn)

    # Print the explanation and the final result as requested.
    print(f"To find a(10), we use symmetry to simplify the problem.")
    print(f"The formula is: a(10) = 4 * (C_straight + 2 * C_turn)")
    print(f"where C_straight is the number of valid 8-step continuations after a 'straight' start ((0,0)->(1,0)->(2,0)),")
    print(f"and C_turn is the number of valid 8-step continuations after a 'turn' start ((0,0)->(1,0)->(1,1)).")
    print("\nCalculating the components:")
    print(f"C_straight = {c_straight}")
    print(f"C_turn = {c_turn}")
    print("\nSubstituting these values into the formula:")
    print(f"a(10) = 4 * ({c_straight} + 2 * {c_turn})")
    print(f"a(10) = 4 * ({c_straight} + {2 * c_turn})")
    print(f"a(10) = 4 * ({c_straight + 2 * c_turn})")
    print(f"a(10) = {total}")

solve_saw_count()
<<<44100>>>