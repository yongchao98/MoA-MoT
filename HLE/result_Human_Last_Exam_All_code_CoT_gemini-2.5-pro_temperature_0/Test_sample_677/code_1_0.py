def find_planar_board_supremum():
    """
    Finds the supremum of the size nm for which the (3,2)-super-knight graph
    on an n x m board (n, m >= 4) is planar.

    The planarity of this graph family is closely related to the inequality:
    3*n*m - 10*n - 10*m + 26 <= 0.
    This script iterates through board dimensions to find the maximum area nm
    that satisfies this condition.
    """
    max_size = 0
    # We only need to check a reasonable range for n and m, as the product
    # will eventually violate the inequality.
    # For m > 10/3, n <= (10m - 26) / (3m - 10). As m -> inf, n -> 10/3.
    # This means n must be small. We can check n up to a reasonable limit,
    # and for each n, find the max m.
    
    print("Searching for the largest rectangle size nm such that 3nm - 10n - 10m + 26 <= 0")
    
    # Based on the inequality n <= (10m - 26)/(3m-10), as m gets large, n approaches 10/3.
    # This means we only need to check small values of n (and by symmetry, m).
    # Let's check n from 4 up to 10, and for each n, find the max m.
    
    best_n, best_m = 0, 0

    for n in range(4, 20): # Check a reasonable range for n
        for m in range(n, 20): # Check m >= n
            # The condition for planarity is 3*n*m - 10*n - 10*m + 26 <= 0
            val = 3 * n * m - 10 * n - 10 * m + 26
            if val <= 0:
                if n * m > max_size:
                    max_size = n * m
                    best_n, best_m = n, m

    print(f"The largest board satisfying the condition is {best_n}x{best_m} with size {max_size}.")
    print(f"Let's verify for ({best_n}, {best_m}):")
    print(f"3 * {best_n} * {best_m} - 10 * {best_n} - 10 * {best_m} + 26 = {3 * best_n * best_m - 10 * best_n - 10 * best_m + 26}")
    
    # Check the next larger sizes to show they fail
    n, m = best_n, best_m + 1
    print(f"Checking board {n}x{m} (size {n*m}):")
    print(f"3 * {n} * {m} - 10 * {n} - 10 * {m} + 26 = {3 * n * m - 10 * n - 10 * m + 26}")

    n, m = 5, 5
    print(f"Checking board {n}x{m} (size {n*m}):")
    print(f"3 * {n} * {m} - 10 * {n} - 10 * {m} + 26 = {3 * n * m - 10 * n - 10 * m + 26}")
    
    print(f"\nThe supremum of the value of nm is {max_size}.")


find_planar_board_supremum()