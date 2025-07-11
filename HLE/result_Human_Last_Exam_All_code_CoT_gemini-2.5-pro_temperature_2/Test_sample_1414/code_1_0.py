import mpmath
import math
from functools import reduce

def find_machin_like_coefficients():
    """
    Solves for the integer coefficients of a given Machin-like formula using the PSLQ algorithm.
    """
    # Set the precision for mpmath. High precision is crucial for PSLQ to succeed.
    mpmath.mp.dps = 150

    # The denominators of the arguments to the arctan functions.
    denominators = [122, 239, 682, 1252, 2855, 12943]

    # Construct the vector of real numbers for which to find an integer relation.
    # The first element is pi/4. The rest are arctan(1/m_i).
    values = [mpmath.pi / 4] + [mpmath.atan(mpmath.mpf(1) / m) for m in denominators]

    # Use the PSLQ algorithm to find a list of integer coefficients `coeffs`
    # such that sum(coeffs[i] * values[i]) is close to zero.
    coeffs = mpmath.pslq(values)

    if coeffs is None:
        print("PSLQ algorithm failed. The precision might be too low.")
        return

    # The found relation is: coeffs[0]*(pi/4) + coeffs[1]*arctan(1/122) + ... = 0
    # The problem's equation is: n*(pi/4) - c1*arctan(1/122) - ... = 0
    # Thus, n = coeffs[0] and c_i = -coeffs[i+1].
    n = coeffs[0]
    c_coeffs = [-k for k in coeffs[1:]]

    # Ensure n is the smallest positive integer.
    # If n from PSLQ is negative, we multiply all coefficients by -1.
    if n < 0:
        n = -n
        c_coeffs = [-k for k in c_coeffs]
    
    # Normalize the coefficients by their greatest common divisor (GCD).
    # PSLQ usually returns a primitive relation, but this is a safeguard.
    all_coeffs_for_gcd = [abs(k) for k in [n] + c_coeffs]
    # In case there's only one non-zero coefficient, gcd would fail with one argument.
    if len(all_coeffs_for_gcd) > 1:
        common_divisor = reduce(math.gcd, all_coeffs_for_gcd)
    else:
        common_divisor = all_coeffs_for_gcd[0] if all_coeffs_for_gcd else 1


    if common_divisor > 1:
        n = n // common_divisor
        c_coeffs = [k // common_divisor for k in c_coeffs]

    c1, c2, c3, c4, c5, c6 = c_coeffs
    
    # --- Output the results ---
    print("The discovered relation is:")
    
    # Build the right-hand side of the equation string
    rhs_terms = []
    denominators_map = {1: 122, 2: 239, 3: 682, 4: 1252, 5: 2855, 6: 12943}
    for i, c_val in enumerate(c_coeffs, 1):
        if c_val != 0:
             rhs_terms.append(f"{c_val}*arctan(1/{denominators_map[i]})")
    
    rhs_string = " + ".join(rhs_terms).replace("+ -", "- ")

    # Print the full equation
    print(f"{n}*pi/4 = {rhs_string}")
    
    # Print the coefficients
    print("\nThe constants are:")
    print(f"n={n}, c1={c1}, c2={c2}, c3={c3}, c4={c4}, c5={c5}, c6={c6}")

if __name__ == '__main__':
    find_machin_like_coefficients()