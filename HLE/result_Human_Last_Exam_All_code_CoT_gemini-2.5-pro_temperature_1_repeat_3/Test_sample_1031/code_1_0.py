def solve_disease_spread():
    """
    Calculates the minimum number of initially infected sites to infect an n x n grid.

    The problem is a case of bootstrap percolation. For an n x n grid with a
    3-neighbor infection threshold, the minimum number of sites required to
    infect the entire grid is given by a known formula:
    - 2n - 1 for n odd
    - 2n - 2 for n even

    This function applies the formula for the given case of n = 14.
    """
    n = 14

    # For n=14, which is an even number, the formula is 2n - 2.
    if n % 2 == 0:
        min_sites = 2 * n - 2
        formula_str = "2 * n - 2"
        calc_str = f"2 * {n} - 2"
    else:
        # This case is not needed for n=14 but included for completeness.
        min_sites = 2 * n - 1
        formula_str = "2 * n - 1"
        calc_str = f"2 * {n} - 1"

    print(f"For an n x n grid where n = {n} (an even number), the minimum number of initial sites is given by the formula: {formula_str}")
    print("The calculation is:")
    print(f"Minimum sites = {calc_str}")
    print(f"Minimum sites = {min_sites}")

solve_disease_spread()