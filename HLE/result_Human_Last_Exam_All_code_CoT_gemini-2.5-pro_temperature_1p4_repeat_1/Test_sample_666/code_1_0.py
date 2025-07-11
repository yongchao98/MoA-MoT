import math

def count_enclosed_poles():
    """
    Counts the number of poles of the function f(z) that lie inside the curve Gamma.
    The curve is described by G(u,v) = 0, where u=x+y, v=x-y.
    A point is inside if G(u,v) < 0.
    The poles are at z = a + 2*pi*k*i, so x=a and y=2*pi*k.
    This gives u = a + 2*pi*k and v = a - 2*pi*k.
    """
    
    def G(u, v):
        """ The equation of the curve in u,v coordinates. """
        try:
            term1 = 3 * (u**2 - 400)**3
            term2 = -20 * u**3 * v**2
            term3 = 1200 * v**4
            return term1 + term2 + term3
        except OverflowError:
            # If numbers are too large, the point is definitely outside.
            return float('inf')

    pole_count = 0
    # The domain for 'a' is given in the problem statement.
    a_domain = range(-2024, 2025)
    
    # We can determine a practical search range for 'k'. The region is bounded.
    # Testing with a=0 shows that for |k|>=3, the point is outside.
    # A safe range for k would be slightly larger.
    k_range = range(-10, 11)

    # We also don't need to check all 'a' values. For large |a|, G(u,v) > 0.
    # A practical range for 'a' can be found by testing, e.g. [-40, 40].
    # But for completeness, we loop through the specified domain of 'a'.
    # We can add an optimization: if for a certain |a|, no k yields an interior
    # point, larger |a| values won't either.
    
    for a in a_domain:
        found_for_this_a = False
        for k in k_range:
            u = a + 2 * math.pi * k
            v = a - 2 * math.pi * k
            
            if G(u, v) < 0:
                pole_count += 1
                found_for_this_a = True

        # Optimization to speed up the counting.
        # If for a given 'a', no 'k' is inside, and 'a' is large enough,
        # we can assume no more poles will be found.
        if not found_for_this_a and abs(a) > 40:
            break
            
    return pole_count

# Calculate the number of poles, N.
N = count_enclosed_poles()

# The integral is 2 * pi * i * N.
integral_coeff = 2 * N

# Output the final equation as requested.
print(f"The number of enclosed poles is N = {N}.")
print(f"The contour integral is given by the equation:oint_Gamma f(z) dz = 2 * pi * i * N")
print(f"Substituting the value of N, we get:oint_Gamma f(z) dz = 2 * pi * i * {N} = {integral_coeff} * pi * i")