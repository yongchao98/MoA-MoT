import sys

# Set a higher recursion limit for deep searches on large boards.
# The default limit (often 1000) is sufficient for N=16, but this is a safeguard.
sys.setrecursionlimit(2000)

class NQueensSolver:
    """
    Solves the N-Queens problem using a backtracking algorithm.
    """
    def __init__(self, n):
        self.n = n
        # board[i] will store the column number of the queen in row i
        self.board = [-1] * n
        self.solution = None

    def is_safe(self, row, col):
        """
        Check if placing a queen at board[row] = col is safe.
        We only need to check for conflicts with queens in previous rows (0 to row-1),
        as queens in subsequent rows have not been placed yet.
        """
        for i in range(row):
            # Check for column conflict
            if self.board[i] == col:
                return False
            # Check for diagonal conflict
            if abs(row - i) == abs(col - self.board[i]):
                return False
        return True

    def solve_util(self, row):
        """
        Utility function to solve the N-Queens problem recursively.
        """
        # Base case: If all queens are placed, we have found a solution
        if row == self.n:
            self.solution = self.board[:]
            return True

        # Try placing a queen in each column of the current row
        for col in range(self.n):
            if self.is_safe(row, col):
                self.board[row] = col
                # Recur to place the rest of the queens
                if self.solve_util(row + 1):
                    return True
                # If placing queen in board[row] = col doesn't lead to a solution,
                # then backtrack (the position is implicitly reset by the next loop iteration).
        
        # If no column works in this row, return False to trigger backtracking
        return False

    def find_solution(self):
        """
        Initializes and runs the solver.
        """
        if self.solve_util(0):
            return self.solution
        else:
            return None

def main():
    """
    Main function to solve the problem for a 16x16 board.
    """
    n = 16
    print(f"Solving for a {n}x{n} chessboard...")
    
    solver = NQueensSolver(n)
    solution_columns = solver.find_solution()

    if solution_columns:
        total_queens = len(solution_columns)
        print(f"\nFound a solution for placing {total_queens} non-attacking queens.")
        
        # We can arbitrarily color half the queens white and half black.
        m = total_queens // 2
        
        print(f"\nThe problem asks for the maximum number m such that m white queens and m black queens")
        print(f"can coexist without attacking each other. This means we are placing 2*m queens in total.")
        print(f"\nThe maximum number of total non-attacking queens on a {n}x{n} board is {total_queens}.")
        print(f"Therefore, the total number of queens, 2*m, must be at most {total_queens}.")
        print(f"This gives the equation: 2 * m = {total_queens}")
        
        print("\nSolving for m:")
        
        num1 = total_queens
        num2 = 2
        result = m
        
        # The final required output format for the equation
        print(f"The equation is: {num1} / {num2} = {result}")

    else:
        print("No solution found for N=16, which contradicts the known theorem.")

if __name__ == "__main__":
    main()