import sys

def solve_n_queens(n):
    """
    Finds a single solution to the N-Queens problem using backtracking.
    Returns a list `p` where p[i] is the column of the queen in row i.
    Returns None if no solution exists.
    """
    # positions[row] = col
    positions = [-1] * n

    def is_safe(row, col):
        """Check if it's safe to place a queen at positions[row] = col."""
        # Check for queen in the same column or diagonals in previous rows
        for r in range(row):
            # Check column conflict
            if positions[r] == col:
                return False
            # Check diagonal conflict
            if abs(row - r) == abs(col - positions[r]):
                return False
        return True

    def find_solution(row):
        """Recursive function to find a solution row by row."""
        # Base case: If all queens are placed, we have a solution
        if row == n:
            return True
        
        # Try placing a queen in each column for the current row
        for col in range(n):
            if is_safe(row, col):
                positions[row] = col
                if find_solution(row + 1):
                    return True
        
        # If no column is safe in this row, backtrack
        return False

    # Start the search from the first row (row 0)
    if find_solution(0):
        return positions
    else:
        return None

def main():
    """
    Solves the problem of finding the maximum number m of coexisting white and black queens.
    """
    N = 16
    
    # Increase recursion limit for the backtracking algorithm, as N=16 is moderately deep.
    sys.setrecursionlimit(N + 5)

    print(f"The problem: Find the maximum number m such that m white queens and m black queens can coexist on a {N}x{N} board without queens of the same color attacking each other.")
    
    print("\nStep 1: Determine the theoretical maximum for m.")
    print(f"The maximum number of non-attacking queens of a single color on an {N}x{N} board is {N}.")
    print(f"Therefore, the value of m cannot be greater than {N}.")

    print(f"\nStep 2: Prove that m = {N} is possible by construction.")
    print("We will find one solution for the 16-Queens problem for the white queens, and then derive a second, disjoint solution for the black queens.")
    
    print("\nFinding a valid placement for 16 white queens...")
    white_queens_pos = solve_n_queens(N)

    if white_queens_pos is None:
        print(f"Error: Could not find a solution for {N}-Queens, which is unexpected.")
        return

    print("Found a solution. This will be the placement for the white queens.")
    
    print("\nStep 3: Create a second, disjoint solution for the black queens.")
    print("We can create a new valid solution by reflecting the white queen positions horizontally.")
    print("If a white queen is at (row, col), the black queen will be at (row, N-1-col).")
    print(f"This is guaranteed to be a disjoint set because the equation 'col = {N-1} - col' (or '2*col = {N-1}') has no integer solution.")

    print("\nStep 4: Display the calculated positions for both sets of queens.")
    print("The following list shows the column for the white queen and the calculated column for the black queen in each row.")

    # Calculate black queen positions
    black_queens_pos = [(N - 1) - c for c in white_queens_pos]
    
    print("-" * 60)
    print(f"{'Row':>4} | {'White Queen Column':>20} | {'Black Queen Column Calculation':>30}")
    print("-" * 60)
    for r in range(N):
        w_col = white_queens_pos[r]
        b_col = black_queens_pos[r]
        equation = f"{b_col:2d} = {N-1} - {w_col:2d}"
        print(f"{r:4d} | {w_col:20d} | {equation:>30}")
    print("-" * 60)

    print("\nConclusion:")
    print(f"We have successfully constructed a valid placement for {N} white queens and {N} black queens.")
    print(f"Since m cannot exceed {N}, the maximum value for m is {N}.")

if __name__ == '__main__':
    main()