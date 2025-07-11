import numpy as np

def solve():
    """
    Calculates the contour integral by finding all poles inside the contour Gamma.
    The integral is 2*pi*i * sum_of_residues.
    The residue at each pole is 1. So we just need to count the poles.
    """

    def is_inside(a, k):
        """
        Checks if the pole z = a + 2*pi*i*k is inside the curve Gamma.
        A point (x,y) is inside if F(x+y, x-y) < 0.
        """
        pi = np.pi
        x = float(a)
        y = 2.0 * pi * k
        u = x + y
        v = x - y

        # From analysis, the curve is confined to 20 <= u <= sqrt(480).
        # We test if the pole's u-coordinate is in this range.
        sqrt_480 = np.sqrt(480.0)
        if not (20.0 <= u <= sqrt_480):
            return False

        # Evaluate the curve's polynomial F(u,v) = 3(u^2-400)^3 - 20u^3v^2 + 1200v^4
        # The interior of the curve corresponds to F(u,v) < 0.
        u2 = u * u
        v2 = v * v
        u2_minus_400 = u2 - 400.0
        
        # We need to handle potential floating point precision issues, but direct evaluation should be sufficient.
        val = 3.0 * (u2_minus_400**3) - 20.0 * (u**3) * v2 + 1200.0 * (v2**2)
        
        return val < 0

    num_poles_inside = 0
    poles_list = []
    
    # The range for 'a' is given in the problem
    a_min = -2024
    a_max = 2024

    # We determine a reasonable search range for k.
    # From 20 <= a + 2*pi*k, we have k >= (20-a)/(2pi).
    # From a + 2*pi*k <= sqrt(480), we have k <= (sqrt(480)-a)/(2pi).
    # For a in [-2024, 2024], k is roughly in [-320, 325].
    # However, analysis shows that v^2 grows much faster than the bounds for large |k|,
    # so we only need to check a small range of k.
    for k in range(-5, 6):
        pi = np.pi
        sqrt_480 = np.sqrt(480.0)
        
        # Determine the range of 'a' for a given 'k' from the condition on u.
        a_lower_bound = np.ceil(20.0 - 2.0 * pi * k)
        a_upper_bound = np.floor(sqrt_480 - 2.0 * pi * k)
        
        a_start = max(a_min, int(a_lower_bound))
        a_end = min(a_max, int(a_upper_bound))

        for a in range(a_start, a_end + 1):
            if is_inside(a, k):
                num_poles_inside += 1
                poles_list.append((a,k))

    print(f"The poles inside the contour are at (a, k): {poles_list}")
    print(f"The total number of poles inside the contour is {num_poles_inside}.")
    
    # The integral is 2 * pi * i * (sum of residues)
    # Sum of residues = number of poles * 1
    sum_of_residues = num_poles_inside
    
    print("\nThe value of the integral is given by the formula:")
    print(f"Integral = 2 * pi * i * (Sum of Residues)")
    print(f"           = 2 * pi * i * {sum_of_residues}")
    print(f"           = {2 * sum_of_residues} * pi * i")

solve()