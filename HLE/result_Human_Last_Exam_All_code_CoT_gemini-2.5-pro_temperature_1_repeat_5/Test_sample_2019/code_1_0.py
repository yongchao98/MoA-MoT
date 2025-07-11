import math

class ToroidalQueensSolver:
    """
    Solves the N-queens problem on a toroidal board for K queens.
    """
    def __init__(self, n, k):
        self.n = n  # Board size
        self.k = k  # Number of queens
        # We can fix the rows due to symmetry. We calculate solutions for one
        # representative row combination, e.g., (0, 1, ..., k-1).
        self.rows = list(range(k))
        self.solutions_for_one_config = 0

    def solve(self):
        """
        Calculates the number of solutions for a fixed set of rows.
        """
        self.backtrack(0, [])
        return self.solutions_for_one_config

    def is_safe(self, queen_idx, new_col, placements):
        """
        Checks if placing a queen at (new_row, new_col) is safe.
        'placements' holds the column positions of queens placed in previous rows.
        """
        new_row = self.rows[queen_idx]

        for i, placed_col in enumerate(placements):
            placed_row = self.rows[i]

            # Check column attack (already handled by design, but good practice)
            if new_col == placed_col:
                return False

            # Check toroidal diagonal attack
            dr = abs(new_row - placed_row)
            dc = abs(new_col - placed_col)
            
            # The toroidal distance is the shorter of the two wrap-around distances
            toroidal_dr = min(dr, self.n - dr)
            toroidal_dc = min(dc, self.n - dc)

            if toroidal_dr == toroidal_dc:
                return False
        return True

    def backtrack(self, queen_idx, placements):
        """
        The backtracking algorithm to find all valid placements.
        """
        # Base case: if all k queens are placed, we found a solution
        if queen_idx == self.k:
            self.solutions_for_one_config += 1
            return

        # Try placing the current queen in each column of the current row
        for col in range(self.n):
            if self.is_safe(queen_idx, col, placements):
                placements.append(col)
                self.backtrack(queen_idx + 1, placements)
                placements.pop()  # Backtrack

def main():
    n = 5  # Size of the board is 5x5
    k = 4  # Number of queens to place

    print(f"To find the number of ways to place {k} non-attacking queens on a {n}x{n} toroidal chessboard:")

    # Step 1: Calculate the number of ways to choose k rows from n
    # This is C(n, k) = n! / (k! * (n-k)!)
    num_row_choices = math.comb(n, k)
    print(f"1. First, we determine the number of ways to choose the {k} rows where the queens will be placed. This is C({n}, {k}) = {num_row_choices}.")

    # Step 2: For one choice of rows, find the number of valid column arrangements
    # Due to symmetry, this number is the same for all choices of rows.
    solver = ToroidalQueensSolver(n, k)
    solutions_per_row_config = solver.solve()
    print(f"2. Second, for any single choice of {k} rows, we find the number of ways to place the queens so they don't attack each other. Using a backtracking search, this number is {solutions_per_row_config}.")
    
    # Step 3: Calculate the total number of ways
    total_ways = num_row_choices * solutions_per_row_config
    print(f"3. Finally, the total number of ways is the product of these two values.")
    print(f"Total ways = (Ways to choose rows) * (Solutions per row choice)")
    print(f"Total ways = {num_row_choices} * {solutions_per_row_config} = {total_ways}")

if __name__ == "__main__":
    main()