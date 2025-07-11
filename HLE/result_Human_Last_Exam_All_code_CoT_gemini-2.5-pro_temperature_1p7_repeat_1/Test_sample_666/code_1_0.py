import numpy as np

def count_poles():
    """
    Counts the number of poles of the function f(z) that lie inside the contour Gamma.
    A pole is at z = a + 2*pi*i*k, corresponding to the point (a, 2*pi*k) in the plane.
    The contour Gamma is contained in the strip |x+y| < sqrt(Y_r),
    where Y_r is the root of K(Y) = -35*Y^3 + 43200*Y^2 - 17280000*Y + 2304000000 = 0.
    The root Y_r is approximately 594.055.
    """
    
    # Root of the polynomial K(Y), found numerically
    Y_r = 594.055333063
    sqrt_Y_r = np.sqrt(Y_r)
    
    a_min_range = -2024
    a_max_range = 2024

    # Determine the range of k
    # The interval [a_min_range + 2*pi*k, a_max_range + 2*pi*k] must overlap with (-sqrt_Y_r, sqrt_Y_r).
    # a_min_range + 2*pi*k < sqrt_Y_r => 2*pi*k < sqrt_Y_r - a_min_range
    # k < (sqrt_Y_r - a_min_range) / (2*pi)
    k_upper_bound = (sqrt_Y_r - a_min_range) / (2 * np.pi) # approx 326.19
    # a_max_range + 2*pi*k > -sqrt_Y_r => 2*pi*k > -sqrt_Y_r - a_max_range
    # k > (-sqrt_Y_r - a_max_range) / (2*pi)
    k_lower_bound = (-sqrt_Y_r - a_max_range) / (2 * np.pi) # approx -326.19
    
    k_min = int(np.ceil(k_lower_bound))
    k_max = int(np.floor(k_upper_bound))

    total_poles = 0
    for k in range(k_min, k_max + 1):
        # For each k, find the range of 'a' that satisfies the condition |a + 2*pi*k| < sqrt_Y_r
        # This is equivalent to -sqrt_Y_r - 2*pi*k < a < sqrt_Y_r - 2*pi*k
        
        lower_a = -sqrt_Y_r - 2 * np.pi * k
        upper_a =  sqrt_Y_r - 2 * np.pi * k
        
        # We need to count integers 'a' in the intersection of three ranges:
        # 1. The function's definition: [-2024, 2024]
        # 2. The geometric constraint from Gamma: (lower_a, upper_a)
        
        # Ceil of lower bound gives the first integer
        start_a = int(np.ceil(lower_a))
        # Floor of upper bound gives the last integer
        end_a = int(np.floor(upper_a))
        
        # Intersect with the given range for 'a'
        final_start_a = max(a_min_range, start_a)
        final_end_a = min(a_max_range, end_a)
        
        if final_end_a >= final_start_a:
            count = final_end_a - final_start_a + 1
            total_poles += count
            
    return total_poles

N = count_poles()
print(f"The number of poles inside the contour is N = {N}")

integral_value_real_part = 0
integral_value_imag_part_coeff = 2 * N

print(f"The integral is 2 * pi * i * N = 2 * pi * i * {N} = {integral_value_imag_part_coeff} * pi * i")
