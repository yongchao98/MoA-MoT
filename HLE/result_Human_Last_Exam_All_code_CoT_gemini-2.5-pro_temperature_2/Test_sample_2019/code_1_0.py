def count_queen_placements():
    """
    Calculates the number of ways to place 4 non-attacking queens
    on a 5x5 toroidal chessboard using a backtracking algorithm.
    """
    # 1. Define problem parameters
    N = 5
    NUM_QUEENS = 4
    
    # We use a list to hold the count, allowing modification within nested functions.
    # A 'nonlocal' keyword or passing it as a parameter would also work.
    solution_count = [0]

    def dist(a, b):
        """Calculates the toroidal distance on a dimension of size N."""
        return min(abs(a - b), N - abs(a - b))

    def is_safe(new_queen, placed_queens):
        """Checks if a new queen at new_queen attacks any of the placed_queens."""
        r2, c2 = new_queen
        for r1, c1 in placed_queens:
            # Check row and column attacks. This is a safeguard; our main
            # strategy for choosing distinct cells implicitly handles this.
            if r1 == r2 or c1 == c2:
                return False
            # Check toroidal diagonal attacks. This is the key constraint.
            if dist(r1, r2) == dist(c1, c2):
                return False
        return True

    def backtrack(k, placements, start_cell):
        """
        Recursive function to find valid placements.
        k: number of queens already placed.
        placements: a list of (row, col) tuples for placed queens.
        start_cell: the linear index of the board cell to start searching from.
                    This enforces a canonical ordering and avoids permutations.
        """
        # Base case: if all queens have been placed, we found a solution.
        if k == NUM_QUEENS:
            solution_count[0] += 1
            return

        # Iterate through possible cells for the next queen, starting from start_cell
        # to ensure we only find unique combinations of cells.
        for i in range(start_cell, N * N):
            r, c = i // N, i % N
            
            # Check if placing a queen at (r, c) is safe
            if is_safe((r, c), placements):
                # If safe, place the queen and recurse for the next one
                placements.append((r, c))
                backtrack(k + 1, placements, i + 1)
                
                # Backtrack: remove the queen to explore other possibilities
                placements.pop()

    # Initial call to start the backtracking process
    # k=0 (placing the first queen), placements=[] (no queens placed yet), start_cell=0 (start from the first square)
    backtrack(0, [], 0)
    
    # Final output as requested
    print("This script finds the number of ways to place Q non-attacking queens on an NxN toroidal board.")
    print(f"Board size (N): {N}")
    print(f"Number of queens (Q): {NUM_QUEENS}")
    print(f"The total number of non-attacking placements is: {solution_count[0]}")


# Run the calculation
count_queen_placements()