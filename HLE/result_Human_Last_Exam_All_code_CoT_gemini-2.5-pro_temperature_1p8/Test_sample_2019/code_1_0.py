def solve_queens_toroidal():
    """
    Calculates the number of ways to place 4 non-attacking queens
    on a 5x5 toroidal chessboard using a backtracking algorithm.
    """
    N = 5
    NUM_QUEENS = 4
    count = 0

    def is_safe(r, c, queens):
        """
        Checks if placing a queen at (r, c) is safe from other queens
        on a toroidal board.
        """
        for qr, qc in queens:
            # Check for row or column conflict
            if r == qr or c == qc:
                return False
            # Check for main diagonal conflict (with wrap-around)
            if (r - c) % N == (qr - qc) % N:
                return False
            # Check for anti-diagonal conflict (with wrap-around)
            if (r + c) % N == (qr + qc) % N:
                return False
        return True

    def backtrack(queens, start_index):
        """
        Recursively explores queen placements to find all solutions.
        
        Args:
            queens (list): A list of (row, col) tuples for placed queens.
            start_index (int): The board square index to start the search from.
                               This ensures each combination is counted only once.
        """
        nonlocal count
        # Base case: if all queens have been placed, we found a solution.
        if len(queens) == NUM_QUEENS:
            count += 1
            return

        # Optimization: stop if not enough squares are left for the remaining queens.
        remaining_squares = N * N - start_index
        if NUM_QUEENS - len(queens) > remaining_squares:
            return

        # Iterate through the board squares from the start_index.
        for i in range(start_index, N * N):
            r = i // N
            c = i % N

            # If the square is safe, place a queen and recurse.
            if is_safe(r, c, queens):
                queens.append((r, c))
                # The next queen search will start from the next square (i + 1).
                backtrack(queens, i + 1)
                # Backtrack: remove the queen to explore other possibilities.
                queens.pop()

    # Start the backtracking process from an empty board at the first square.
    backtrack([], 0)
    print(count)

if __name__ == "__main__":
    solve_queens_toroidal()