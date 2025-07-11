def solve_disease_spread():
    """
    Calculates the minimum number of initially infected sites to infect a 14x14 grid.

    The method is based on a modified checkerboard pattern.
    1. A checkerboard pattern where (i+j) is even guarantees all 2x2 blocks are hit.
       This is a necessary condition to prevent stable healthy blocks.
       The number of such sites is n*n / 2.
    2. However, a pure checkerboard pattern leaves two corners with an odd coordinate sum healthy.
       These corners only have 2 neighbors, so they can never become infected.
    3. To fix this, we must also infect these two corners from the start.
    4. The total number of sites is the sum of the checkerboard sites and the two extra corners.
    """
    n = 14
    # Number of sites where i+j is even
    checkerboard_sites = (n * n) // 2
    # The corners (0, 13) and (13, 0) have odd coordinate sums (13).
    # The other two corners (0,0) and (13,13) have even sums and are already included.
    # These two corners must be added to the initial set to ensure they become infected.
    extra_corners = 2
    
    total_sites = checkerboard_sites + extra_corners
    
    print(f"For an {n}x{n} grid:")
    print(f"Number of sites in the base checkerboard pattern ('i+j' is even): {checkerboard_sites}")
    print(f"Number of additional corner sites needed: {extra_corners}")
    print(f"The minimum number of sites is the sum of these two values.")
    print(f"Final Equation: {checkerboard_sites} + {extra_corners} = {total_sites}")
    print(f"The minimum number of initially infected sites is {total_sites}.")

solve_disease_spread()