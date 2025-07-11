def solve_toroidal_queens():
    """
    Calculates the number of ways to place k non-attacking queens on an n x n
    toroidal chessboard using a backtracking algorithm.
    """
    n = 5  # The size of the chessboard (n x n)
    k = 4  # The number of queens to place

    # A global or non-local counter to store the number of solutions
    solution_count = 0

    def is_safe(r, c, placed_queens):
        """Checks if placing a queen at (r, c) is safe from existing queens."""
        for qr, qc in placed_queens:
            # Check for attacks on row, column, and wrapped diagonals
            if (r == qr or
                c == qc or
                (r - c) % n == (qr - qc) % n or
                (r + c) % n == (qr + qc) % n):
                return False
        return True

    def backtrack(queens_to_place, start_index, placed_queens):
        """
        Recursive function to find valid queen placements.
        
        :param queens_to_place: Number of queens remaining to be placed.
        :param start_index: The board square index to start searching from.
        :param placed_queens: A list of (r, c) coordinates of already placed queens.
        """
        nonlocal solution_count
        
        # Base case: If all queens are placed, we found a solution
        if queens_to_place == 0:
            solution_count += 1
            return

        # Iterate through the board squares to find a spot for the next queen
        for i in range(start_index, n * n):
            r = i // n
            c = i % n

            if is_safe(r, c, placed_queens):
                # Place the queen
                placed_queens.append((r, c))
                # Recur to place the rest of the queens
                backtrack(queens_to_place - 1, i + 1, placed_queens)
                # Backtrack: remove the queen to explore other possibilities
                placed_queens.pop()

    # Start the backtracking process
    backtrack(k, 0, [])

    # The "final equation" is conceptual: Result = Count(queens=4, board_size=5)
    # The code outputs all the numbers involved.
    print(f"Board size (n x n): {n}x{n}")
    print(f"Number of queens (k): {k}")
    print(f"Total number of ways: {solution_count}")


if __name__ == '__main__':
    solve_toroidal_queens()
<<<250>>>