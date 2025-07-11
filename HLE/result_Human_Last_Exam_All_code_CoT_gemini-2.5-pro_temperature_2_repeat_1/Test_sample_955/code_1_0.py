import collections

def solve_grid_problem():
    """
    Calculates the value of n for the given grid problem.
    The function iterates through even values of n, calculates the number of
    reachable cells for each n using a Breadth-First Search (BFS), and
    checks if the probability condition (66%) is met.
    """

    def count_reachable_cells(n):
        """
        Calculates the number of unique cells reachable within max_moves.

        Args:
            n: The size of the n x n grid.

        Returns:
            The total number of reachable cells.
        """
        # The grid must be large enough for the starting cell 'c2' (3, 2)
        if n < 3:
            return 0

        start_pos = (3, 2)
        max_moves = 3

        # BFS queue stores tuples of (position, moves_taken)
        queue = collections.deque([(start_pos, 0)])

        # A dictionary to store the minimum moves to reach each cell,
        # which also acts as a 'visited' set.
        distances = {start_pos: 0}

        def is_valid(x, y):
            return 1 <= x <= n and 1 <= y <= n

        while queue:
            current_pos, current_moves = queue.popleft()

            if current_moves >= max_moves:
                continue

            x, y = current_pos
            new_moves = current_moves + 1

            # --- Move Type 1: Diagonal Moves ---
            # For each of the four diagonal directions
            for dx, dy in [(1, 1), (1, -1), (-1, 1), (-1, -1)]:
                nx, ny = x, y
                # Trace along the entire diagonal path
                while True:
                    nx += dx
                    ny += dy
                    if not is_valid(nx, ny):
                        break  # Stop at the grid edge

                    new_pos = (nx, ny)
                    # Add to queue if this is a new cell or a shorter path is found
                    if new_pos not in distances or distances[new_pos] > new_moves:
                        distances[new_pos] = new_moves
                        queue.append((new_pos, new_moves))

            # --- Move Type 2: Border Moves ---
            is_on_border = (x == 1 or x == n or y == 1 or y == n)
            if is_on_border:
                # Check all four adjacent neighbors
                for dx, dy in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
                    ax, ay = x + dx, y + dy

                    if not is_valid(ax, ay):
                        continue
                    
                    # The adjacent cell must also be on the border
                    is_adj_on_border = (ax == 1 or ax == n or ay == 1 or ay == n)
                    if is_adj_on_border:
                        new_pos = (ax, ay)
                        if new_pos not in distances or distances[new_pos] > new_moves:
                            distances[new_pos] = new_moves
                            queue.append((new_pos, new_moves))
        
        return len(distances)

    # Main loop to find the correct value of n
    # Start with n=4, the smallest even integer that can contain 'c2'
    n = 4
    while True:
        num_reachable = count_reachable_cells(n)
        total_cells = n * n

        # Check for the condition: num_reachable / total_cells == 0.66
        # To avoid floating point issues, use the equivalent integer equation:
        # 50 * num_reachable == 33 * total_cells
        if 50 * num_reachable == 33 * total_cells:
            print(f"Solution found for n = {n}")
            print(f"Total grid cells: {total_cells}")
            print(f"Cells reachable within 3 moves: {num_reachable}")
            # Print the equation with all the numbers
            print(f"Final Equation: {num_reachable} / {total_cells} = {num_reachable/total_cells}")
            print(f"The value of n is {n}.")
            return n

        # Break if n gets unreasonably large to prevent infinite loops
        if n > 100:
            print("Could not find a solution for n up to 100.")
            return None
        
        n += 2  # n must be even

# Run the solver
solve_grid_problem()