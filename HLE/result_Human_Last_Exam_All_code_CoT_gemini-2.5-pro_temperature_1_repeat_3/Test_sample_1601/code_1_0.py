import numpy as np
from scipy.integrate import quad

def calculate_omega_measure():
    """
    Calculates the area of the set Omega based on an analytical approximation
    of the separatrix between blowing-up and non-blowing-up solutions.
    """

    # The domain for initial conditions (a0, b0) is [-10, 1] x [10, 20].
    # The approximate condition for the desired blow-up is a0 * (a0^2 + a0 - b0^2 / 2) > 0.

    # Case 1: a0 is in (0, 1]
    # We need a0^2 + a0 > b0^2 / 2.
    # Max(LHS) = 1^2 + 1 = 2. Min(RHS) = 10^2 / 2 = 50.
    # The condition 2 > 50 is never met.
    area_a0_pos = 0.0

    # Case 2: a0 is in [-1, 0)
    # We need a0^2 + a0 - b0^2 / 2 < 0, since a0 < 0.
    # LHS = a0(a0+1) is negative. RHS = b0^2/2 is positive.
    # The condition is always met.
    # Area = width * height = (0 - (-1)) * (20 - 10)
    area_a0_neg1_to_0 = 1.0 * 10.0

    # Case 3: a0 is in [-10, -1)
    # We need a0^2 + a0 - b0^2 / 2 < 0, which means b0 > sqrt(2 * (a0^2 + a0)).
    # Let f(a0) = sqrt(2 * (a0^2 + a0)).
    
    # Find intersection of the boundary curve f(a0) with b0 = 10.
    # a0^2 + a0 - 50 = 0
    a_intersect = (-1 - np.sqrt(1**2 - 4*1*(-50))) / 2.0
    
    # The area calculation is split into two parts at a_intersect.
    # Part A: a0 from a_intersect to -1.
    # Here f(a0) < 10, so the condition b0 > f(a0) is true for all b0 in [10, 20].
    # Area = width * height = ((-1) - a_intersect) * (20 - 10)
    area_a0_large_neg_A = (-1.0 - a_intersect) * 10.0
    
    # Part B: a0 from -10 to a_intersect.
    # Here f(a0) > 10. The blow-up region is where b0 is between f(a0) and 20.
    # Area is given by the integral of (20 - f(a0)) da0.
    integrand = lambda a: 20.0 - np.sqrt(2 * (a**2 + a))
    area_a0_large_neg_B, _ = quad(integrand, -10.0, a_intersect)
    
    # Total area is the sum of these parts.
    total_area = area_a0_pos + area_a0_neg1_to_0 + area_a0_large_neg_A + area_a0_large_neg_B

    print("The estimated size of the set Omega, m(Omega), is calculated by summing the areas from different subdomains of a0:")
    print(f"1. For a0 in (0, 1]: The blow-up condition is not met. Area = {area_a0_pos:.2f}")
    print(f"2. For a0 in [-1, 0): The condition is always met. Area = {area_a0_neg1_to_0:.2f}")
    print(f"3. For a0 in [-10, -1): The area is split at a0 = {a_intersect:.2f}")
    print(f"   - For a0 in [{a_intersect:.2f}, -1): Area = {area_a0_large_neg_A:.2f}")
    print(f"   - For a0 in [-10, {a_intersect:.2f}): Area from integral = {area_a0_large_neg_B:.2f}")
    
    print("\nFinal Equation for m(Omega):")
    print(f"{area_a0_pos:.2f} + {area_a0_neg1_to_0:.2f} + {area_a0_large_neg_A:.2f} + {area_a0_large_neg_B:.2f} = {total_area:.2f}")
    
calculate_omega_measure()