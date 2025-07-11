def solve_infection_problem():
    """
    Calculates the minimum number of initially infected sites to infect a 14x14 grid.

    The problem is equivalent to finding the minimum size of a hitting set for all
    2x2 squares on the grid.
    """
    n = 14

    # The logic relies on the finding that for an n x n grid where n is even,
    # the minimum number of sites to hit all 2x2 squares is (n/2)^2.
    # This is established by a lower bound from a disjoint tiling and an
    # upper bound from a specific construction.

    # Calculation for n = 14
    side_count = n / 2
    min_sites = side_count * side_count

    print("The grid size is n x n, where n = 14.")
    print("The problem is equivalent to finding the minimum number of sites to hit every 2x2 square.")
    print("For an even-sized grid, this number can be calculated as (n/2)^2.")
    print("The final calculation is:")
    print(f"({n} / 2) * ({n} / 2) = {int(side_count)} * {int(side_count)} = {int(min_sites)}")

solve_infection_problem()