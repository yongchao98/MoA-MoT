import math
from functools import reduce
import mpmath

def gcd_list(numbers):
    """Computes the greatest common divisor of a list of integers."""
    # Take absolute values and filter out any zeros
    num_list = [abs(int(x)) for x in numbers if x != 0]
    if not num_list:
        return 1
    return reduce(math.gcd, num_list)

def solve_arctan_relation():
    """
    Finds the integer coefficients for the given Machin-like formula.
    """
    # Set a high precision for the calculations to ensure correctness
    mpmath.mp.dps = 200

    # The denominators of the arctangent arguments
    y_values = [122, 239, 682, 1252, 2855, 12943]

    # The equation to solve is:
    # c1*atan(1/y1) + ... + c6*atan(1/y6) - n*(pi/4) = 0
    # We construct a vector for PSLQ corresponding to this relation.
    terms = []
    for y in y_values:
        terms.append(mpmath.atan(mpmath.mpf(1) / y))
    terms.append(-mpmath.pi / 4)

    # Run the PSLQ algorithm to find the integer relation vector
    # The resulting coefficients will correspond to [c1, c2, ..., c6, n]
    coeffs = mpmath.pslq(terms, tol=1e-180)
    
    # Convert coefficients to integers
    coeffs = [int(round(c)) for c in coeffs]

    # PSLQ finds a relation, but it might not be primitive.
    # Divide by the GCD to get the smallest integer coefficients.
    common_divisor = gcd_list(coeffs)
    if common_divisor > 1:
        coeffs = [c // common_divisor for c in coeffs]

    # Extract the coefficients c_i and n
    c = coeffs[:-1]
    n = coeffs[-1]

    # We need the smallest positive n. If n is negative, we flip the signs
    # of all coefficients, which corresponds to multiplying the equation by -1.
    if n < 0:
        n = -n
        c = [-val for val in c]

    # Print the full equation with the found coefficients
    print("The final equation is:")
    equation_str = f"{n} * pi/4 = "
    first_term = True
    for i in range(len(c)):
        if c[i] == 0:
            continue
        
        term_val = c[i]
        term_str = f"{abs(term_val)}*arctan(1/{y_values[i]})"
        
        if not first_term:
            equation_str += " + " if term_val > 0 else " - "
        elif term_val < 0:
            equation_str += "-"

        equation_str += term_str
        first_term = False
    print(equation_str)

    # Print the coefficients individually
    print("\nThe solution for the constants is:")
    print(f"n = {n}")
    for i in range(len(c)):
        print(f"c_{i+1} = {c[i]}")

    # Return the final answer as a comma-separated string for submission
    final_answer_string = f"{n},{','.join(map(str, c))}"
    return final_answer_string

if __name__ == '__main__':
    solve_arctan_relation()
