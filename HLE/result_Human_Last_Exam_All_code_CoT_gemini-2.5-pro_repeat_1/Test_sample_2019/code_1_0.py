def solve_queens():
    """
    Calculates the number of ways to place 4 non-attacking queens on a 5x5 toroidal board.
    """
    N = 5  # Size of the board
    Q = 4  # Number of queens

    def is_safe(pos, placements):
        """
        Checks if placing a queen at `pos` is safe with respect to `placements`.
        A position is safe if it's not on the same row, column, or toroidal diagonal
        as any queen already in the placements list.
        """
        r_new, c_new = pos
        for r_old, c_old in placements:
            # Check row and column
            if r_new == r_old or c_new == c_old:
                return False
            # Check toroidal diagonals using modulo arithmetic
            if (r_new - c_new) % N == (r_old - c_old) % N:
                return False
            if (r_new + c_new) % N == (r_old + c_old) % N:
                return False
        return True

    def count_solutions(queens_to_place, start_index, placements):
        """
        Recursively finds the number of ways to place the remaining queens.
        - queens_to_place: How many more queens we need to place.
        - start_index: The linear index (0-24) to start searching from. This prevents
                       duplicate solutions by ensuring a fixed placement order.
        - placements: A list of (row, col) tuples for queens already placed.
        """
        # Base case: If all queens have been placed, we found one valid solution.
        if queens_to_place == 0:
            return 1

        count = 0
        # Iterate through all cells starting from the specified index.
        for i in range(start_index, N * N):
            r = i // N
            c = i % N

            # If it's safe to place a queen here...
            if is_safe((r, c), placements):
                # ...place it and recurse to find solutions for the remaining queens.
                # The next search will start from the next cell (i + 1).
                new_placements = placements + [(r, c)]
                count += count_solutions(queens_to_place - 1, i + 1, new_placements)
        
        return count

    # Start the search: place 4 queens, starting from cell 0, with no queens placed yet.
    total_ways = count_solutions(Q, 0, [])
    print(total_ways)

solve_queens()
<<<250>>>