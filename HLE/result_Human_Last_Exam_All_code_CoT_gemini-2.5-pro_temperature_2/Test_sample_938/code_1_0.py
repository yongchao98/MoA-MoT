import numpy as np
import math

def solve():
    """
    Solves the problem by finding the roots of the characteristic polynomial
    and calculating the limit based on their moduli.
    """
    # The characteristic equation is r^4 - 88r - 57 = 0.
    # The coefficients are [1, 0, 0, -88, -57] for r^4, r^3, r^2, r^1, r^0.
    coeffs = [1, 0, 0, -88, -57]
    
    # Find the roots of the polynomial.
    roots = np.roots(coeffs)
    
    # We expect one large positive real root (alpha), one real root between -1 and 0,
    # and a pair of complex conjugate roots.
    # The modulus of a complex number z = a + bj is sqrt(a^2 + b^2).
    # numpy.abs() calculates the modulus for both real and complex numbers.
    moduli = np.abs(roots)
    
    # Sort the moduli in descending order to identify alpha and rho.
    sorted_moduli = sorted(moduli, reverse=True)
    
    alpha = sorted_moduli[0]
    rho = sorted_moduli[1]
    
    # The problem asks for the limit L = lim_{n->inf} ln(s_n)/n.
    # My analysis showed that the leading asymptotic term of s_n is of the order (alpha*rho)^n.
    # Therefore, L = ln(alpha*rho) = ln(alpha) + ln(rho).
    
    limit_val = math.log(alpha) + math.log(rho)
    
    # We need to find the integer part of 10^4 * limit_val.
    result = 10000 * limit_val
    
    integer_part = int(result)
    
    # The problem asks to output the equation to be calculated
    # Let's show the numerical values for clarity
    print(f"The characteristic equation is r^4 - 88*r - 57 = 0.")
    # Printing roots can be complicated due to formatting, so we show the moduli.
    print(f"The moduli of the roots, sorted, are: {sorted_moduli[0]:.6f}, {sorted_moduli[1]:.6f}, {sorted_moduli[2]:.6f}, {sorted_moduli[3]:.6f}")
    print(f"The dominant root modulus is alpha = {alpha:.6f}")
    print(f"The second largest root modulus is rho = {rho:.6f}")
    print(f"The value to calculate is the integer part of 10000 * (ln(alpha) + ln(rho)).")
    print(f"10000 * (ln({alpha:.6f}) + ln({rho:.6f})) = 10000 * ({math.log(alpha):.6f} + {math.log(rho):.6f}) = {result:.6f}")
    print(f"The integer part is {integer_part}.")
    
    
solve()
