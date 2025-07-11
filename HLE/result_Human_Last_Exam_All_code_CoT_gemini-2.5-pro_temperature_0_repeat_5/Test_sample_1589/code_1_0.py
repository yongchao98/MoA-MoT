import sys

# It's recommended to increase the recursion limit for deep searches like N-Queens.
sys.setrecursionlimit(2000)

class Solver:
    """
    This class solves the problem by finding the number of valid permutations for a given n.
    """
    def __init__(self, n):
        self.n = n
        self.count = 0
        # sigma[i] will store the column of the queen in row i
        self.sigma = [0] * (n + 1)
        # Helper arrays for the N-Queens backtracking search
        self.used_cols = [False] * (n + 1)
        # For diagonal checks: i - sigma[i] and i + sigma[i] must be unique
        self.used_diag1 = [False] * (2 * n)
        self.used_diag2 = [False] * (2 * n + 1)

    def check_additional_condition(self):
        """
        Checks condition (c): sigma(j) - sigma(i) != i + j for all i, j.
        This is equivalent to abs(sigma(j) - sigma(i)) != i + j for i < j,
        because sigma(i) - sigma(j) = i + j is impossible as LHS is negative and RHS is positive.
        We only need to check for pairs (i, j) where a violation is possible,
        which is when i + j <= n - 1.
        """
        for i in range(1, self.n + 1):
            for j in range(i + 1, self.n + 1):
                # A violation can only occur if i+j is in the range of |sigma(j)-sigma(i)|,
                # which is [1, n-1]. So we check for i+j <= n-1.
                if i + j <= self.n - 1:
                    if abs(self.sigma[j] - self.sigma[i]) == i + j:
                        return False
        return True

    def search(self, row):
        """
        Backtracking search to find N-Queens solutions.
        """
        if row > self.n:
            # A solution to the N-Queens problem is found (sigma).
            # Now, check the additional problem-specific condition.
            if self.check_additional_condition():
                self.count += 1
            return

        for col in range(1, self.n + 1):
            # Check if placing a queen at (row, col) is safe (no attacks)
            # diag1_idx: row - col, shifted to be non-negative
            # diag2_idx: row + col
            diag1_idx = row - col + self.n - 1
            diag2_idx = row + col
            if (not self.used_cols[col] and
                    not self.used_diag1[diag1_idx] and
                    not self.used_diag2[diag2_idx]):
                
                self.sigma[row] = col
                self.used_cols[col] = True
                self.used_diag1[diag1_idx] = True
                self.used_diag2[diag2_idx] = True

                self.search(row + 1)

                # Backtrack: undo the choice for the current row
                self.used_cols[col] = False
                self.used_diag1[diag1_idx] = False
                self.used_diag2[diag2_idx] = False

    def solve(self):
        """
        Starts the search and returns the final count.
        """
        self.search(1)
        return self.count

if __name__ == '__main__':
    n = 10
    solver = Solver(n)
    result = solver.solve()
    # The problem asks for a(10), which is the final calculated number.
    # The final equation is a(10) = result.
    print(f"a({n}) = {result}")
