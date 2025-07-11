import mpmath

def solve_machin_like_formula():
    """
    Finds the integer coefficients for a Machin-like formula for pi/4
    using the PSLQ algorithm.
    """
    # Set the precision for the calculations.
    # The PSLQ algorithm requires high precision to find the correct integer relations.
    mpmath.mp.dps = 100

    # The given denominators for the arctan arguments.
    x_values = [122, 239, 682, 1252, 2855, 12943]

    # The problem is to find an integer relation of the form:
    # n*(pi/4) = c1*arctan(1/x1) + c2*arctan(1/x2) + ... + c6*arctan(1/x6)
    # This can be rewritten as a homogeneous linear combination equal to zero:
    # c1*arctan(1/x1) + ... + c6*arctan(1/x6) - n*(pi/4) = 0
    # We can use the PSLQ algorithm to find the integer coefficients [c1, ..., c6, -n].

    # We create a vector of the high-precision real numbers.
    vec = [mpmath.atan(mpmath.mpf(1)/x) for x in x_values]
    vec.append(mpmath.pi / 4)

    # Run the PSLQ algorithm to find the integer relation.
    coeffs = mpmath.pslq(vec)

    # The result of pslq is a list of integers [c1, c2, ..., c6, -n].
    # We extract the coefficients from the result.
    c = coeffs[:-1]
    n = -coeffs[-1]

    # The problem asks for the smallest positive n.
    # If PSLQ returns a relation with a negative n, we can multiply all
    # coefficients by -1 to make n positive, and the equation remains valid.
    if n < 0:
        n = -n
        c = [-ci for ci in c]

    # Print the found coefficients.
    print(f"Found coefficients for the smallest positive n:")
    print(f"n = {n}")
    for i in range(len(c)):
        print(f"c_{i+1} = {c[i]}")

    # Print the full equation with the found coefficients for verification.
    # The prompt requires to output each number in the final equation.
    print("\nThe equation is:")
    equation_parts = []
    for i in range(len(c)):
        equation_parts.append(f"({c[i]})*arctan(1/{x_values[i]})")
    equation_str = " + ".join(equation_parts)
    print(f"{n} * pi/4 = {equation_str}")
    
    # Provide the final answer in the required format.
    final_answer_string = f"{n},{','.join(map(str, c))}"
    print(f"\n<<<" + final_answer_string + ">>>")

solve_machin_like_formula()