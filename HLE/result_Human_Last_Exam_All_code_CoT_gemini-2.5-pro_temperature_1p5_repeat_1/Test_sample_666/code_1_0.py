import math

def solve():
    """
    This script calculates the contour integral by counting the number of poles
    of the function f(z) that are enclosed by the curve Gamma.
    """
    
    pole_count = 0
    pi = math.pi
    
    # Based on the analysis of the curve's equation, we know that any
    # enclosed pole a + 2*pi*i*k must satisfy |a| < 20 and |k| <= 3.
    # So we iterate through these ranges.
    for a in range(-19, 20):
        for k in range(-3, 4):
            x = float(a)
            y = 2.0 * pi * k
            
            u = x + y
            w = x - y
            
            # A pole at (x, y) is inside the curve Gamma if the following inequality holds.
            # This inequality is derived from the corrected and simplified equation of the curve.
            # ((x+y)^2 - 400)^3 + 400(x-y)^4 < 0
            if (u**2 - 400.0)**3 + 400.0 * w**4 < 0:
                pole_count += 1
                
    # According to the Residue Theorem, the integral is 2 * pi * i * (sum of residues).
    # Since the residue at each enclosed pole is 1, the sum is equal to the pole count.
    integral_value_coefficient = 2 * pole_count
    
    print("The number of enclosed poles is calculated to be:", pole_count)
    print("The value of the contour integral is given by the equation:")
    print(f"oint_Gamma f(z) dz = 2 * pi * i * (Number of enclosed poles)")
    print(f"oint_Gamma f(z) dz = 2 * pi * i * {pole_count}")
    print(f"oint_Gamma f(z) dz = {integral_value_coefficient} * pi * i")

solve()