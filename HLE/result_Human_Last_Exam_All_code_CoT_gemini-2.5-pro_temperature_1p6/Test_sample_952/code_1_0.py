def solve_diamond_problem():
    """
    Calculates the largest number of diamonds that can be placed on a 2024x2024 grid
    such that for any valid arrangement, at least one diamond can be moved.

    The logic is based on finding the size of the smallest "frozen" configuration,
    where no diamond can be moved. A frozen configuration can be made by placing
    diamonds on all cells of a single color on a checkerboard, which for an
    N x N grid (with N even) is (N*N)/2.

    The largest value `k` for which *every* configuration is movable must be
    one less than the size of this frozen configuration.
    """
    N = 2024
    
    # Total number of cells
    total_cells = N * N
    
    # Number of diamonds in a frozen checkerboard configuration
    max_diamonds_frozen = total_cells // 2
    
    # The largest value k must be one less than this
    largest_k = max_diamonds_frozen - 1
    
    print(f"The size of the table is N = {N}.")
    print(f"A 'frozen' arrangement can be made with (N * N) / 2 diamonds.")
    print(f"The calculation is: ({N} * {N}) / 2 - 1 = {largest_k}")
    print(f"The largest value is {largest_k}.")

solve_diamond_problem()