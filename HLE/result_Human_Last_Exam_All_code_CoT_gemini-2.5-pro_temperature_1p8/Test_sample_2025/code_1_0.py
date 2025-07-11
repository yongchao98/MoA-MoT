import numpy as np
from scipy.integrate import quad

def get_conditional_pdf(z, x):
    """
    Calculates the conditional PDF f_{Z|X1}(z|x).
    This function implements the analytical solution to the integral:
    integral from d to M of (1/y) * 6 * F(y) * (1-F(y)) * f(y) dy
    where d = |z-x|, M = max(x, 1-x), and F(y), f(y) are the CDF and PDF
    of the distance |X-x| for X ~ U[0,1].
    """
    d = np.abs(z - x)
    m = np.min([x, 1 - x])
    M = np.max([x, 1 - x])
    
    # This term appears in both cases
    xx = x * (1 - x) # m - m^2 or (1-m)- (1-m)^2, depends on m=x or m=1-x

    # If d is outside the feasible range of distances, density is 0
    if d > M:
        return 0.0

    # Case 1: The integration range [d, M] is entirely in the region where y > m
    if d >= m:
        # We integrate 6*(-y + (1-2m) + (xx)/y) dy from d to M
        val_at_M = -M**2/2 + (1-2*m)*M + xx*np.log(M) if M > 0 else 0
        val_at_d = -d**2/2 + (1-2*m)*d + xx*np.log(d) if d > 0 else -np.inf
        # handle log(0) case for d=0, the limit is 0 if xx=0 (x=0 or x=1), else diverges. 
        # But d>=m and d=0 implies m=0, so x=0 or x=1.
        if d==0 and xx > 0: # only for z=x, with x in (0,1). d=0 <= m. Handled below.
           val_at_d = 0 # this is a limit, let's treat it.
        elif d==0 and xx==0:
           val_at_d = 0
        if np.isinf(val_at_d): # a failsafe for d=0, m>0 which is handled below but to be safe
             val_at_d=0 # it should be covered by d<m

        result = 6 * (val_at_M - val_at_d)
        
    # Case 2: The integration range [d, M] is split into [d, m] and [m, M]
    else: # d < m
        # Part 1: Integral from d to m
        # We integrate 24*(1-2y) dy
        part1 = 24 * (m - m**2 - (d - d**2))
        
        # Part 2: Integral from m to M
        val_at_M = -M**2/2 + (1-2*m)*M + xx*np.log(M) if M > 0 else 0
        val_at_m = -m**2/2 + (1-2*m)*m + xx*np.log(m) if m > 0 else 0
        
        # log(0) case for m=0 (x=0 or x=1)
        if m==0 and xx>0: # Should not happen
            val_at_m = 0
        elif m==0 and xx==0: # This happens at x=0, x=1. m=0. xx=0. 
            val_at_m = 0

        part2 = 6 * (val_at_M - val_at_m)
        result = part1 + part2
        
    return result

def calculate_final_pdf(z_val):
    """
    Calculates the final PDF f_Z(z_val) by integrating the conditional PDF.
    Splits the integral into regions for numerical stability and correctness.
    """
    # For z=0.2, the condition d=|0.2-x| <= m=min(x, 1-x) holds for x in [0.1, 0.6]
    # The condition d > m holds for x in [0, 0.1) U (0.6, 1]
    
    # Integrand for the numerical integration
    integrand = lambda x: get_conditional_pdf(z_val, x)

    # Calculate the integral over each region
    # Region 1: [0, 0.1]
    integral1, _ = quad(integrand, 0, 0.1)
    
    # Region 2: [0.1, 0.6]
    integral2, _ = quad(integrand, 0.1, 0.6)

    # Region 3: [0.6, 1.0]
    integral3, _ = quad(integrand, 0.6, 1.0)
    
    # The total value is the sum of the integrals
    total_value = integral1 + integral2 + integral3
    
    print(f"The PDF f_Z(z) is the integral of a conditional PDF f(z|x) over x from 0 to 1.")
    print(f"For z = {z_val}, we split the integral into three parts:")
    print(f"Integral from 0.0 to 0.1: {integral1:.6f}")
    print(f"Integral from 0.1 to 0.6: {integral2:.6f}")
    print(f"Integral from 0.6 to 1.0: {integral3:.6f}")
    print(f"\nSumming these values gives the final result:")
    print(f"f_Z({z_val}) = {integral1:.6f} + {integral2:.6f} + {integral3:.6f} = {total_value:.6f}")
    
    return total_value

if __name__ == '__main__':
    z_value = 0.2
    final_result = calculate_final_pdf(z_value)

