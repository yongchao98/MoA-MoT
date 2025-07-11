import mpmath
import math
from functools import reduce

def find_gcd_of_list(numbers):
    """
    Computes the greatest common divisor of a list of integers.
    """
    # Filter out zeros and take absolute values.
    non_zero_numbers = [abs(x) for x in numbers if x != 0]
    if not non_zero_numbers:
        # This case should not be reached in this problem.
        return 1
    return reduce(math.gcd, non_zero_numbers)

def solve_machin_like_formula():
    """
    Finds the integer coefficients for the given Machin-like formula using the PSLQ algorithm.
    """
    # Set the precision for the calculation. PSLQ needs high precision.
    mpmath.mp.dps = 100

    # The denominators of the arctan arguments from the problem statement
    denominators = [122, 239, 682, 1252, 2855, 12943]

    # The vector of values for which to find an integer relation.
    # The relation is of the form: k_0*v_0 + k_1*v_1 + ... = 0
    # Our equation is n*pi/4 - c_1*arctan(1/y1) - ... = 0
    # So v_0 = pi/4, v_1 = arctan(1/y1), ...
    # and we will find k_0 = n, k_1 = -c_1, ...
    print("Constructing the vector of values...")
    vector = [mpmath.pi / 4] + [mpmath.atan(mpmath.mpf(1) / d) for d in denominators]

    # Use the PSLQ algorithm to find the integer coefficients [k_0, k_1, ..., k_6]
    print("Running PSLQ algorithm to find integer relation...")
    coefficients = mpmath.pslq(vector)

    # Reduce the coefficients by their greatest common divisor (GCD)
    # to get the smallest integer relation.
    common_divisor = find_gcd_of_list(coefficients)
    reduced_coeffs = [c // common_divisor for c in coefficients]

    # Extract n and c_i from the result based on the equation structure
    n = reduced_coeffs[0]
    cs = [-c for c in reduced_coeffs[1:]]

    # The problem asks for the smallest positive n.
    # If n is negative, we can multiply the entire relation by -1.
    if n < 0:
        n = -n
        cs = [-c for c in cs]

    # Unpack the coefficients for printing
    c1, c2, c3, c4, c5, c6 = cs

    # Build and print the equation string to show the result
    print("\nFound the following relation:")
    equation_parts = [f"{n}*pi/4 ="]
    first_term = True
    for c, d in zip(cs, denominators):
        if c == 0:
            continue
        
        sign = ""
        # For subsequent terms, the sign is always explicit
        if not first_term:
            sign = "+ " if c > 0 else "- "
        # For the first term, only show sign if it's negative
        elif c < 0:
            sign = "- "
        
        abs_c = abs(c)
        
        if abs_c == 1:
            term = f"arctan(1/{d})"
        else:
            term = f"{abs_c}*arctan(1/{d})"
            
        equation_parts.append(f"{sign}{term}")
        first_term = False

    # The strip() removes potential trailing space
    print(" ".join(equation_parts))

    # Print the final coefficients
    print("\nThe integer coefficients are:")
    print(f"n = {n}")
    print(f"c1 = {c1}")
    print(f"c2 = {c2}")
    print(f"c3 = {c3}")
    print(f"c4 = {c4}")
    print(f"c5 = {c5}")
    print(f"c6 = {c6}")
    
    # Print the answer in the specified format for parsing
    final_answer = f"{n},{c1},{c2},{c3},{c4},{c5},{c6}"
    print(f"\nFinal answer in the format n,c1,c2,c3,c4,c5,c6:\n{final_answer}")


if __name__ == '__main__':
    solve_machin_like_formula()
