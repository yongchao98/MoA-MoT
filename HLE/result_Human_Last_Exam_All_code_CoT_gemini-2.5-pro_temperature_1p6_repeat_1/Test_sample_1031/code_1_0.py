import math

def solve_disease_spread(n):
    """
    Calculates the minimum number of initially infected sites to infect an n x n grid.
    The rule is a healthy site becomes infected if it has at least 3 infected neighbors.

    Args:
        n (int): The dimension of the grid.

    Returns:
        int: The minimum number of initial sites.
    """
    if n <= 0:
        return 0

    # The base strategy is a checkerboard pattern.
    # We infect all sites (i, j) where i + j is even.
    num_checkerboard_sites = math.ceil(n * n / 2.0)
    
    # This pattern infects all non-corner sites effectively.
    # However, corner sites can only be infected if they are in the initial set.
    # We must check if the checkerboard pattern includes all four corners.
    
    num_missing_corners = 0
    
    # Corner (0, 0). 0 + 0 = 0 (even). Included in the checkerboard pattern.
    
    # For n even, like n=14:
    # Corner (0, n-1): 0 + (n-1) is odd. Missing.
    # Corner (n-1, 0): (n-1) + 0 is odd. Missing.
    # Corner (n-1, n-1): (n-1) + (n-1) is even. Included.
    if n % 2 == 0:
        num_missing_corners = 2
        # For n even, the number of black/white squares is exactly n*n/2
        num_checkerboard_sites = n * n / 2
        
    # For n odd, all corners have an even sum of coordinates and are included.
    # e.g., n=15: (0,0), (0,14), (14,0), (14,14) all have even sums.

    min_sites = int(num_checkerboard_sites + num_missing_corners)
    
    print(f"For an {n}x{n} grid:")
    print(f"Number of sites in the base checkerboard pattern: {int(num_checkerboard_sites)}")
    print(f"Number of corners missed by the pattern: {num_missing_corners}")
    print("The final calculation is:")
    print(f"{int(num_checkerboard_sites)} + {num_missing_corners} = {min_sites}")

# The specific problem is for n = 14
n = 14
solve_disease_spread(n)
