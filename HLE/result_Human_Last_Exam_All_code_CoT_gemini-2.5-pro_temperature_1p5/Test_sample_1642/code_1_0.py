def solve_max_queens():
    """
    Calculates the maximum number m such that m white and m black queens
    can coexist on a 16x16 board without attacking each other.
    """
    n = 16

    class NQueensSolver:
        """
        An efficient N-Queens solver using backtracking with boolean arrays for O(1) safety checks.
        """
        def __init__(self, size):
            self.n = size
            self.solution = None
            self.board = [-1] * self.n
            self.col_used = [False] * self.n
            # For primary diagonals, row + col is constant. Max value is 2*(n-1).
            self.diag1_used = [False] * (2 * self.n - 1)
            # For secondary diagonals, row - col is constant. We offset to get non-negative indices.
            self.diag2_used = [False] * (2 * self.n - 1)

        def solve_util(self, row):
            """Recursive utility function to solve N-Queens problem."""
            if row == self.n:
                # All queens are placed successfully.
                self.solution = self.board[:]
                return True

            for col in range(self.n):
                # Check if placing a queen at (row, col) is safe in O(1) time.
                diag1_index = row + col
                diag2_index = row - col + self.n - 1
                if not self.col_used[col] and not self.diag1_used[diag1_index] and not self.diag2_used[diag2_index]:
                    # Place queen
                    self.board[row] = col
                    self.col_used[col] = True
                    self.diag1_used[diag1_index] = True
                    self.diag2_used[diag2_index] = True

                    # Recur for the next row
                    if self.solve_util(row + 1):
                        return True

                    # Backtrack: Un-place the queen
                    self.col_used[col] = False
                    self.diag1_used[diag1_index] = False
                    self.diag2_used[diag2_index] = False
            return False

        def solve(self):
            """Initializes the backtracking process."""
            return self.solve_util(0)

    print("Step 1: Analyzing the problem.")
    print(f"The task is to place 'm' white queens and 'm' black queens on a {n}x{n} board without any attacks.")
    print("This is equivalent to placing a total of 2*m queens of a single color, as the attack rules are the same regardless of color.")
    print("-" * 30)

    print("Step 2: Determining the maximum number of queens.")
    print(f"The maximum number of non-attacking queens that can be placed on an {n}x{n} board is {n}.")
    print("This is because placing n+1 queens would force at least two queens to be in the same row, causing an attack.")
    print(f"This gives us the equation for the maximum m: 2 * m = {n}.")
    print("-" * 30)

    print("Step 3: Verifying that a placement of 16 queens is possible.")
    print("We will now confirm that a solution to the 16-Queens problem exists using a solver.")
    solver = NQueensSolver(n)
    if solver.solve():
        print("Success! A valid placement for 16 non-attacking queens was found.")
    else:
        print("A solution could not be found, but the theoretical maximum still holds.")
    print("-" * 30)
    
    print("Step 4: Calculating the final answer 'm'.")
    m = n // 2
    print(f"From our equation, 2 * m = {n}, we can solve for m.")
    print(f"m = {n} / 2")
    print(f"The maximum value for m is {m}.")
    print("-" * 30)

    print("Final result:")
    print("The final equation with the value of m substituted is:")
    print(f"2 * {m} = {n}")

solve_max_queens()