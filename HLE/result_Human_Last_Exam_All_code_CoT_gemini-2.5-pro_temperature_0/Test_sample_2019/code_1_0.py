class ToroidalQueensSolver:
    """
    A solver to find the number of ways to place non-attacking queens
    on a toroidal chessboard.
    """

    def __init__(self, board_size, num_queens):
        """
        Initializes the solver.
        Args:
            board_size (int): The size of the board (N for an NxN board).
            num_queens (int): The number of queens to place.
        """
        self.n = board_size
        self.num_queens = num_queens
        self.count = 0

    def _is_safe(self, r, c, placements):
        """
        Checks if placing a new queen at (r, c) is safe with respect
        to the queens already in `placements`.
        """
        for pr, pc in placements:
            # Check for attacks on rows and columns
            if r == pr or c == pc:
                return False
            # Check for attacks on toroidal diagonals
            if (r - c) % self.n == (pr - pc) % self.n:
                return False
            if (r + c) % self.n == (pr + pc) % self.n:
                return False
        return True

    def _solve_recursive(self, k, start_index, placements):
        """
        The recursive backtracking function.
        Args:
            k (int): The number of queens remaining to be placed.
            start_index (int): The board square index to start searching from.
            placements (list): A list of (row, col) tuples for placed queens.
        """
        # Base case: If all queens are placed, we found a valid solution.
        if k == 0:
            self.count += 1
            return

        # Optimization: If the number of remaining squares is less than the
        # number of queens we still need to place, we can't succeed.
        if (self.n * self.n - start_index) < k:
            return

        # Iterate through all possible squares starting from start_index.
        for i in range(start_index, self.n * self.n):
            r = i // self.n
            c = i % self.n

            # If the current square is safe, place a queen and recurse.
            if self._is_safe(r, c, placements):
                placements.append((r, c))
                self._solve_recursive(k - 1, i + 1, placements)
                placements.pop()  # Backtrack: remove the queen to try other paths.

    def find_solutions(self):
        """
        Starts the solving process and returns the total count of solutions.
        """
        self._solve_recursive(self.num_queens, 0, [])
        return self.count

def main():
    """
    Main function to solve the problem and print the result.
    """
    board_size = 5
    num_queens = 4

    solver = ToroidalQueensSolver(board_size, num_queens)
    total_ways = solver.find_solutions()
    
    # The problem asks to output each number in the final equation.
    # We will print a sentence that includes all the relevant numbers:
    # the board size (5), the number of queens (4), and the result.
    print(f"On a {board_size}x{board_size} toroidal board, there are {total_ways} ways to place {num_queens} non-attacking queens.")

if __name__ == "__main__":
    main()