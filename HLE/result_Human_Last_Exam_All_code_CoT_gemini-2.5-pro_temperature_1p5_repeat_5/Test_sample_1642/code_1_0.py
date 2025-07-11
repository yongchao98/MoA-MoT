import sys

def solve_m_queens_problem():
    """
    This function determines and explains the solution to the problem of placing
    m white and m black queens on an NxN board.
    """
    N = 16
    
    # Set a higher recursion limit for larger N like 16.
    # The default limit is often 1000, which might be insufficient.
    sys.setrecursionlimit(N + 5)

    print(f"Determining the maximum number m for an {N}x{N} board.")
    print("-" * 60)
    
    # 1. Interpretation
    print("Step 1: Understanding the Problem")
    print("The goal is to find the maximum `m` where `m` white queens and `m` black queens can coexist.")
    print("This means we need to find two sets of queen positions, one for white and one for black, such that:")
    print("  - Within each set, no two queens attack each other.")
    print("  - The two sets of positions are disjoint (no square is shared).")
    print("-" * 60)
    
    # 2. Upper Bound
    print("Step 2: Finding the Maximum Possible Value for m")
    print(f"For a single set of queens on an {N}x{N} board, the maximum number of non-attacking queens is {N}.")
    print("Placing more than N would require at least two queens to share a row, which is not allowed.")
    print(f"Therefore, the theoretical maximum for m is {N}.")
    print("-" * 60)
    
    # 3. Test m = 16
    print(f"Step 3: Testing if m = {N} is Possible")
    print(f"To achieve m = {N}, we must find two disjoint solutions to the {N}-Queens problem.")
    print("Let's first find one solution for the white queens using a backtracking algorithm.")

    # A solution is represented by an array `p` where p[row] = column.
    white_queen_cols = [-1] * N
    
    def is_safe(board, row, col):
        """Checks if placing a queen at board[row]=col is safe from previous queens."""
        for i in range(row):
            # Check for column conflict and diagonal conflict
            if board[i] == col or abs(board[i] - col) == abs(i - row):
                return False
        return True

    def find_one_solution_util(board, row):
        """Recursive utility to find one solution."""
        if row == N:
            return True  # All queens are placed successfully
        
        for col in range(N):
            if is_safe(board, row, col):
                board[row] = col
                if find_one_solution_util(board, row + 1):
                    return True
        return False

    if not find_one_solution_util(white_queen_cols, 0):
        print("Error: Could not find a solution for 16-Queens. This should not happen.")
        return

    print("Successfully found a base solution for the first set of queens (e.g., White Queens).")
    print("-" * 60)

    # 4. Construct the second solution and prove disjointness
    print("Step 4: Constructing and Verifying a Second, Disjoint Solution")
    print("We can create a second valid solution by horizontally reflecting the first solution.")
    print(f"If the white queen in row `i` is in column `c`, the black queen in row `i` will be in column `({N-1} - c)`.")
    
    black_queen_cols = [(N - 1) - c for c in white_queen_cols]
    
    print("\nThe positions for the White Queens are (row, col):")
    # Using format specifiers in an f-string to output the final equation parts.
    print(' '.join(f"({r},{c})" for r, c in enumerate(white_queen_cols)))
    
    print("\nThe positions for the Black Queens (reflected) are (row, col):")
    print(' '.join(f"({r},{c})" for r, c in enumerate(black_queen_cols)))

    print("\n\nTo be disjoint, a white queen and a black queen can never be on the same square.")
    print("This means for any row `i`, `white_queen_cols[i]` must not equal `black_queen_cols[i]`.")
    print("The condition for them to be on the same square would be: `c = (N-1) - c`.")
    print(f"This simplifies to the equation: 2 * c = {N-1}")
    print(f"For our {N}x{N} board, this is: 2 * c = 15")
    print("Since `c` (the column index) must be an integer, there is no solution to this equation.")
    print("Therefore, the original solution and its reflection are always disjoint.")
    print("-" * 60)

    # 5. Conclusion
    print("Conclusion:")
    print(f"We have successfully constructed two disjoint sets of {N} non-attacking queens.")
    print(f"This proves that m = {N} is possible.")
    print(f"Since m cannot be greater than {N}, the maximum value of m is {N}.")

if __name__ == '__main__':
    solve_m_queens_problem()