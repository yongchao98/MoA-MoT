import math

def solve_integral():
    """
    Calculates the contour integral by counting the poles of f(z) inside the contour Gamma.

    The method assumes the contour Gamma is the circle x^2 + y^2 = 200, based on
    identifying four points on the curve that lie on this circle.

    The integral is 2 * pi * i * N, where N is the number of poles inside Gamma.
    """
    R_squared = 200
    a_min, a_max_range = -2024, 2024
    
    # Determine the maximum |k| to check
    k_max = int(math.sqrt(R_squared / (4 * math.pi**2)))
    
    total_poles = 0
    pole_counts_by_k = {}
    
    # Count poles for k=0
    a_sq_max_k0 = R_squared
    a_max_k0 = int(math.sqrt(a_sq_max_k0))
    # Ensure 'a' is within its defined range [-2024, 2024]
    if a_max_k0 > a_max_range:
        a_max_k0 = a_max_range
    
    count_k0 = 2 * a_max_k0 + 1
    total_poles += count_k0
    pole_counts_by_k[0] = count_k0
    
    # Count poles for k > 0 and k < 0
    for k in range(1, k_max + 1):
        a_sq_max = R_squared - (2 * math.pi * k)**2
        if a_sq_max < 0:
            continue
        
        a_max_val = int(math.sqrt(a_sq_max))
        # Ensure 'a' is within its defined range
        if a_max_val > a_max_range:
            a_max_val = a_max_range

        # Count for +k and -k
        count_for_k = 2 * (2 * a_max_val + 1)
        total_poles += count_for_k
        pole_counts_by_k[k] = 2 * a_max_val + 1
        pole_counts_by_k[-k] = 2 * a_max_val + 1
        
    print(f"The number of poles inside the contour for k=0 is: {pole_counts_by_k[0]}")
    for k in range(1, k_max + 1):
        print(f"The number of poles inside the contour for k=+/
-{k} is: {pole_counts_by_k[k]} for each")
    
    # The final equation for the number of poles
    pole_sum_str_parts = [str(pole_counts_by_k[0])]
    for k in range(1, k_max + 1):
        pole_sum_str_parts.append(str(2 * pole_counts_by_k[k]))

    print(f"The total number of poles is the sum: {' + '.join(pole_sum_str_parts)} = {total_poles}")

    # The integral value is 2 * pi * i * N
    integral_val = 2j * math.pi * total_poles
    
    print("\nThe value of the contour integral is N * 2 * pi * i:")
    print(f"Final calculation: {total_poles} * 2 * {math.pi} * 1j = {integral_val}")

solve_integral()
