class ToroidalQueensSolver:
    """
    A class to solve the N-queens problem on a toroidal board.
    """

    def __init__(self, board_size, num_queens):
        """
        Initializes the solver.
        
        Args:
            board_size (int): The size of the square board (N for an NxN board).
            num_queens (int): The number of queens to place.
        """
        self.board_size = board_size
        self.num_queens = num_queens
        self.solution_count = 0

    def is_safe(self, r, c, placed_queens):
        """
        Checks if placing a queen at (r, c) is safe from already placed queens.
        
        Args:
            r (int): The row of the new queen.
            c (int): The column of the new queen.
            placed_queens (list): A list of (row, col) tuples for existing queens.
            
        Returns:
            bool: True if the placement is safe, False otherwise.
        """
        for qr, qc in placed_queens:
            # Check for row or column attacks
            if r == qr or c == qc:
                return False
            
            # Check for toroidal diagonal attacks
            row_dist = abs(r - qr)
            col_dist = abs(c - qc)
            
            # The toroidal distance is the minimum of the direct and wrap-around distances
            toroidal_row_dist = min(row_dist, self.board_size - row_dist)
            toroidal_col_dist = min(col_dist, self.board_size - col_dist)
            
            if toroidal_row_dist == toroidal_col_dist:
                return False
        return True

    def solve(self, queens_to_place, start_index, placed_queens):
        """
        Recursively finds valid placements using backtracking.
        
        Args:
            queens_to_place (int): Number of queens remaining to be placed.
            start_index (int): The linear index of the square to start searching from.
                               This prevents generating permutations of the same placement.
            placed_queens (list): List of (row, col) tuples for placed queens.
        """
        # Base case: If all queens are placed, we found a solution.
        if queens_to_place == 0:
            self.solution_count += 1
            return

        # Iterate through all board squares, starting from start_index.
        for i in range(start_index, self.board_size * self.board_size):
            r = i // self.board_size
            c = i % self.board_size
            
            # Check if it's safe to place a queen here
            if self.is_safe(r, c, placed_queens):
                # Place the queen
                placed_queens.append((r, c))
                # Recurse to place the next queen, starting from the next square
                self.solve(queens_to_place - 1, i + 1, placed_queens)
                # Backtrack: remove the queen to explore other possibilities
                placed_queens.pop()

    def find_solutions(self):
        """
        Starts the backtracking search and prints the final result.
        """
        self.solve(self.num_queens, 0, [])
        # The final "equation" is the statement of the result.
        # The numbers in the equation are the board size and the number of queens.
        print(f"The number of ways to place {self.num_queens} non-attacking queens on a {self.board_size}x{self.board_size} toroidal chessboard is: {self.solution_count}")


if __name__ == '__main__':
    # Define problem parameters
    BOARD_SIZE = 5
    NUM_QUEENS = 4
    
    # Create a solver and run it
    solver = ToroidalQueensSolver(BOARD_SIZE, NUM_QUEENS)
    solver.find_solutions()