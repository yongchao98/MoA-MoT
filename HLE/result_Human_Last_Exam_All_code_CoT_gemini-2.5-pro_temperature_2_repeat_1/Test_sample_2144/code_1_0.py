import numpy as np
import cmath

def solve_trajectory():
    """
    Solves for the position x0 where y(x0) = -3.
    """
    # The problem reduces to solving the quartic equation p^4 - 18p - 27 = 0 for p.
    # The coefficients of the polynomial are [1, 0, 0, -18, -27].
    coeffs = [1, 0, 0, -18, -27]
    roots = np.roots(coeffs)

    print("The equation to solve for parameter p is: p^4 - 18p - 27 = 0")
    print(f"The roots of this equation are: {roots}")

    # The domain for p requires 2p + 3 > 0, which means p > -1.5.
    # We need to find real roots that satisfy this condition.
    valid_p_values = []
    for p in roots:
        if abs(p.imag) < 1e-9 and p.real > -1.5:
            valid_p_values.append(p.real)
    
    print(f"\nValid real roots for p (where p > -1.5) are: {valid_p_values}")

    # For each valid p, calculate the corresponding x coordinate.
    # The formula for x is x(p) = -(3p + 3) / sqrt(2p + 3)
    final_answer = None
    for p_val in valid_p_values:
        # Check if p is close to an integer, likely the intended solution
        if abs(p_val - round(p_val)) < 1e-9:
            p_int = round(p_val)
            print(f"\nFound an integer root p = {p_int}. Calculating x0 for this value.")
            
            numerator = -(3 * p_int + 3)
            denominator_val = 2 * p_int + 3
            denominator = cmath.sqrt(denominator_val)
            x0 = numerator / denominator
            
            print("Calculation of x0:")
            print(f"p = {p_int}")
            print(f"x0 = -(3*p + 3) / sqrt(2*p + 3)")
            print(f"x0 = -({3*p_int} + {3}) / sqrt({2*p_int} + {3})")
            print(f"x0 = {numerator} / sqrt({denominator_val})")
            print(f"x0 = {numerator} / {denominator.real}")
            print(f"x0 = {x0.real}")
            final_answer = x0.real
            break
            
    # In case the integer root is not found first, just process the first valid root found
    if final_answer is None and valid_p_values:
        print("\nNo integer root found, presenting the first valid solution.")
        p_val = valid_p_values[0] # Take the first one found
        numerator = -(3 * p_val + 3)
        denominator_val = 2 * p_val + 3
        denominator = cmath.sqrt(denominator_val)
        x0 = numerator / denominator
        final_answer = x0.real

solve_trajectory()

print("<<<",-4.0,">>>")