import math

def solve():
    """
    Calculates the contour integral by finding the number of poles of f(z)
    inside the curve Gamma, which is determined to be the circle |z|^2 = 200.
    """
    R_squared = 200
    a_min = -2024
    a_max = 2024
    
    num_poles = 0
    
    # Case k = 0
    # a^2 < 200
    max_a_val_sq = R_squared
    max_a = math.sqrt(max_a_val_sq)
    # Number of integers a such that |a| < max_a is 2*floor(max_a) + 1
    # since max_a is not an integer.
    count_a = 2 * math.floor(max_a) + 1
    num_poles += count_a
    
    # Case k > 0
    k = 1
    while True:
        k_term = (2 * math.pi * k)**2
        if k_term >= R_squared:
            break
        
        max_a_val_sq = R_squared - k_term
        max_a = math.sqrt(max_a_val_sq)
        
        # Count for a in [-floor(max_a), floor(max_a)]
        # This check is actually not needed as max_a < 15, well within [-2024, 2024]
        count_a_k = 2 * math.floor(max_a) + 1
        
        # Add counts for both +k and -k
        num_poles += 2 * count_a_k
        k += 1
        
    integral_coeff = 2 * num_poles
    
    print(f"The curve Gamma is identified as the circle |z|^2 = 200.")
    print(f"The number of poles of f(z) inside Gamma is N = {num_poles}.")
    print(f"The integral is 2 * pi * i * N = 2 * {num_poles} * pi * i = {integral_coeff} * pi * i.")
    
solve()