import math

def count_poles_in_circle(R_squared):
    """
    Counts the number of poles of the form z = a + 2*pi*i*k inside
    the circle |z|^2 < R_squared.
    The integer 'a' is restricted to the range [-2024, 2024], but this
    range is wide enough not to be a limiting factor.
    """
    pole_count = 0
    pi = math.pi
    
    # Case k = 0
    # a^2 < R_squared
    max_a_k0 = math.floor(math.sqrt(R_squared - 1e-9)) # small epsilon for stability
    pole_count += (2 * max_a_k0 + 1)
    
    # Case k > 0 (will be multiplied by 2 for k < 0)
    k = 1
    while True:
        k_term = (2 * pi * k)**2
        if k_term >= R_squared:
            break
        
        max_a_squared = R_squared - k_term
        max_a = math.floor(math.sqrt(max_a_squared))
        
        # Number of 'a' values for a given k > 0 is 2*max_a + 1
        # We multiply by 2 to account for both k and -k
        pole_count += 2 * (2 * max_a + 1)
        k += 1
        
    return pole_count

def solve_integral():
    """
    Calculates the contour integral based on the analysis.
    """
    print("Step 1: The residue of f(z) at each pole is 1.")
    print("Step 2: The integral is 2 * pi * i * N, where N is the number of poles inside Gamma.")
    print("Step 3: Based on testing points, the curve Gamma is assumed to enclose the same poles as the circle |z|^2 = 200.")

    R_squared = 200
    N = count_poles_in_circle(R_squared)
    
    print(f"\nStep 4: Counting the poles inside |z|^2 < {R_squared}.")
    print(f"The number of poles inside the contour is N = {N}")
    
    integral_coeff = 2 * N
    
    print("\nFinal equation for the integral:")
    print(f"Integral = N * 2 * pi * i")
    print(f"Integral = {N} * 2 * pi * i")
    print(f"Integral = {integral_coeff} * pi * i")

solve_integral()