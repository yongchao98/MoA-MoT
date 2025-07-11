class ToroidalNQueens:
    """
    Calculates the number of ways to place K non-attacking queens on an N x N toroidal chessboard.
    """

    def __init__(self, n, k):
        """
        Initializes the solver.
        :param n: The size of the board (n x n).
        :param k: The number of queens to place.
        """
        self.n = n
        self.k = k
        self.count = 0
        self.board_size = n * n

    def _is_safe(self, r1, c1, queens):
        """
        Checks if placing a queen at (r1, c1) is safe with respect to already placed queens.
        A queen at (r1, c1) attacks (r2, c2) if:
        - r1 = r2 (same row)
        - c1 = c2 (same column)
        - (r1 - c1) % n = (r2 - c2) % n (same diagonal)
        - (r1 + c1) % n = (r2 + c2) % n (same anti-diagonal)
        """
        for r2, c2 in queens:
            if r1 == r2 or c1 == c2:
                return False
            if (r1 - c1) % self.n == (r2 - c2) % self.n:
                return False
            if (r1 + c1) % self.n == (r2 + c2) % self.n:
                return False
        return True

    def _solve_recursive(self, queens_to_place, start_index, queens):
        """
        The recursive backtracking function to find valid placements.
        """
        # Base case: if all queens are placed, we found a valid solution.
        if queens_to_place == 0:
            self.count += 1
            return

        # Optimization: if remaining squares are not enough to place remaining queens
        if self.board_size - start_index < queens_to_place:
            return

        # Iterate through all squares on the board starting from 'start_index'
        for i in range(start_index, self.board_size):
            r = i // self.n
            c = i % self.n

            # If it's safe to place a queen here
            if self._is_safe(r, c, queens):
                # Place the queen
                queens.append((r, c))
                # Recur for the next queen, starting from the next square
                self._solve_recursive(queens_to_place - 1, i + 1, queens)
                # Backtrack: remove the queen to explore other possibilities
                queens.pop()

    def find_solutions(self):
        """
        Public method to start the search and return the total count.
        """
        self.count = 0
        if self.k == 0:
            return 1
        if self.k > self.board_size:
            return 0
        
        # Start the recursive search
        self._solve_recursive(self.k, 0, [])
        return self.count

# --- Main Execution ---
# Define the problem parameters
N_board_size = 5
K_queens = 4

# Create a solver instance and run the calculation
solver = ToroidalNQueens(n=N_board_size, k=K_queens)
result = solver.find_solutions()

# Print the final result including the numbers from the problem statement
print(f"On a {N_board_size}x{N_board_size} toroidal board, {K_queens} non-attacking queens can be placed in {result} ways.")