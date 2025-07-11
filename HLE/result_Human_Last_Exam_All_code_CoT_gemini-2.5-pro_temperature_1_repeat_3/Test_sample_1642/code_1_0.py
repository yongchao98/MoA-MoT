import sys

# It's possible the recursion depth could be an issue for N=16, so increase the limit.
sys.setrecursionlimit(2000)

class QueenSolver:
    """
    A solver to find one solution for the N-Queens problem.
    """
    def __init__(self, n):
        self.n = n
        self.solution = None

    def is_safe(self, board, row, col):
        """
        Checks if it's safe to place a queen at board[row][col].
        This function is called when "col" queens are already placed in
        columns from 0 to col-1. So we need to check only left side for
        attacking queens.
        """
        for i in range(row):
            # Check for same column or diagonals
            if board[i] == col or abs(board[i] - col) == abs(i - row):
                return False
        return True

    def solve_util(self, board, row):
        """A recursive utility function to solve N Queens problem."""
        # Base case: If all queens are placed then return true
        if row >= self.n:
            self.solution = board[:] # Make a copy of the solution
            return True

        # Consider this column and try placing this queen in all rows one by one
        for col in range(self.n):
            if self.is_safe(board, row, col):
                board[row] = col
                # Recur to place rest of the queens
                if self.solve_util(board, row + 1):
                    return True
                # If placing queen in board[i][col] doesn't lead to a solution,
                # then remove queen from board[i][col] (backtrack)
                # No explicit removal needed as we just overwrite board[row] in the next iteration.
        return False

    def find_solution(self):
        """Finds one solution and returns it."""
        board = [-1] * self.n
        if not self.solve_util(board, 0):
            return None
        # Convert to 1-based indexing for clarity
        return [col + 1 for col in self.solution]

def main():
    """
    Main function to solve the problem and print the output.
    """
    N = 16
    print("Thinking Process:")
    print("1. The problem asks for the maximum m where m white and m black queens can coexist.")
    print("2. This is interpreted as finding two non-attacking queen configurations of size m that are disjoint.")
    print(f"3. The maximum possible value for m on a {N}x{N} board is {N}.")
    print(f"4. We can prove m={N} is achievable by finding two disjoint {N}-Queens solutions.")
    print("-" * 20)
    print("Finding a solution for the first set of queens (White)...")
    
    solver = QueenSolver(N)
    white_queen_solution = solver.find_solution()

    if white_queen_solution is None:
        print(f"Could not find a solution for {N}-Queens.")
        return

    print("Found a solution for the White Queens.")
    print("The positions (row, column) are:")
    white_coords = []
    for row_idx, col_idx in enumerate(white_queen_solution):
        white_coords.append(f"({row_idx + 1}, {col_idx})")
    print(", ".join(white_coords))
    print("-" * 20)

    print("Constructing a second, disjoint solution for the Black Queens...")
    print("The construction rule for the column of the black queen in row 'r' is:")
    print(f"black_col(r) = ({N} + 1) - white_col(r)")
    
    black_queen_solution = [(N + 1) - col for col in white_queen_solution]
    
    print("The positions (row, column) for the Black Queens are:")
    black_coords = []
    for row_idx, col_idx in enumerate(black_queen_solution):
        black_coords.append(f"({row_idx + 1}, {col_idx})")
    print(", ".join(black_coords))
    print("-" * 20)

    print("Conclusion:")
    print(f"We have successfully placed {N} white queens and {N} black queens on the board")
    print("in two valid, non-overlapping configurations.")
    print(f"Therefore, the maximum number m is {N}.")
    print("\nFinal Equation:")
    print(f"m = {N}")

if __name__ == "__main__":
    main()
