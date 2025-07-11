import math

def calculate_min_max_cells(R):
    """
    Calculates the minimal and maximal number of grid cells a circle of radius R can cross.
    """
    R2 = R * R

    # Calculate r_2(R^2) for N_max
    # R^2 = 500^2 = 250000 = 2^4 * 5^6
    # The odd prime factor is 5 (form 4k+1) with power 6. No 4k+3 factors.
    # r_2(R^2) = 4 * (6 + 1) = 28
    r2_R2 = 28

    # N_max occurs when the center is near a grid vertex (e.g., (eps, eps))
    # Formula: N_max = 8R - 2*r_2(R^2) + 4
    max_cells = 8 * R - 2 * r2_R2 + 4

    # Calculate K for N_min
    # N_min occurs when the center is in the middle of a cell (e.g., (0.5, 0.5))
    # K = sum_{k=1 to R} is_square(floor(sqrt(R^2 - k^2)))
    K = 0
    for k in range(1, R + 1):
        val = math.sqrt(R2 - k * k)
        floor_val = math.floor(val)
        
        # Check if floor_val is a non-negative integer that is a perfect square
        if floor_val >= 0:
            sqrt_floor_val = math.isqrt(floor_val)
            if sqrt_floor_val * sqrt_floor_val == floor_val:
                K += 1
    
    # Formula: N_min = 8R - 4K
    min_cells = 8 * R - 4 * K

    print(f"For a circle of radius R = {R}:")
    print(f"The number of integer points on the circle x^2 + y^2 = R^2 is r_2(R^2) = {r2_R2}.")
    print(f"The maximal number of cells crossed is 8 * {R} - 2 * {r2_R2} + 4 = {max_cells}.")
    print(f"The auxiliary value K is {K}.")
    print(f"The minimal number of cells crossed is 8 * {R} - 4 * {K} = {min_cells}.")
    print(f"The answer is ({min_cells}, {max_cells}).")

# Run the calculation for R=500
calculate_min_max_cells(500)