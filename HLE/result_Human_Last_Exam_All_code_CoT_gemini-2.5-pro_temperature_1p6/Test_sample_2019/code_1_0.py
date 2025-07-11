import math
import itertools

def solve_queens_problem():
    """
    Calculates the number of ways to place 4 non-attacking queens on a 5x5
    toroidal chessboard using a combinatorial method.
    """
    N = 5  # Board size
    K = 4  # Number of queens

    # Step 1: Calculate the number of ways to choose K rows and K columns.
    # We use the math.comb function, which is equivalent to n! / (k! * (n-k)!).
    ways_to_choose_rows = math.comb(N, K)
    ways_to_choose_cols = math.comb(N, K)

    # Step 2: For a fixed KxK subgrid (e.g., rows/cols 0..K-1), count the
    # number of valid queen placements (permutations) that satisfy the
    # toroidal non-attacking conditions on the NxN board.
    
    # We use a representative set of rows and columns for the KxK subgrid.
    subgrid_indices = range(K)
    
    valid_perms_count = 0
    # Iterate through all permutations of column assignments for the rows.
    for p in itertools.permutations(subgrid_indices):
        # A permutation 'p' maps row `i` to column `p[i]`.
        # Queens are at (0, p[0]), (1, p[1]), ..., (K-1, p[K-1]).
        
        # Check diagonal attack conditions using modulo N arithmetic.
        # The values of (row - col) % N must be unique for all queens.
        # The values of (row + col) % N must also be unique for all queens.
        diffs = {(i - p[i]) % N for i in subgrid_indices}
        sums = {(i + p[i]) % N for i in subgrid_indices}

        # If the number of unique differences and sums equals K, the permutation is valid.
        if len(diffs) == K and len(sums) == K:
            valid_perms_count += 1

    # Step 3: The total number of ways is the product of the three parts.
    total_ways = ways_to_choose_rows * ways_to_choose_cols * valid_perms_count

    # Print the final equation with each number, as requested.
    print(f"The total number of ways is the product of:")
    print(f"- Ways to choose {K} rows from {N}: C({N},{K}) = {ways_to_choose_rows}")
    print(f"- Ways to choose {K} columns from {N}: C({N},{K}) = {ways_to_choose_cols}")
    print(f"- Valid non-attacking placements on a {K}x{K} subgrid: {valid_perms_count}")
    print(f"\nFinal Equation: {ways_to_choose_rows} * {ways_to_choose_cols} * {valid_perms_count} = {total_ways}")


solve_queens_problem()
<<<50>>>