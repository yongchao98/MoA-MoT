import math

def solve_equation():
    """
    Calculates and prints the four roots of the given polynomial in increasing order.
    The roots were determined analytically by factoring the polynomial.
    """
    # The four roots of the polynomial are sqrt(14), 2*sqrt(6), sqrt(34), and 2*sqrt(11).
    # We will represent them symbolically and calculate their numerical values.
    roots_data = [
        ("sqrt(14)", math.sqrt(14)),
        ("2*sqrt(6)", 2 * math.sqrt(6)),
        ("sqrt(34)", math.sqrt(34)),
        ("2*sqrt(11)", 2 * math.sqrt(11))
    ]

    # Sort the list based on the numerical value (the second element of the tuple)
    sorted_roots = sorted(roots_data, key=lambda x: x[1])

    # The problem asks to output each root of the equation.
    print("The four roots of the equation, in increasing order, are:")
    for symbolic_form, numeric_value in sorted_roots:
        print(f"{symbolic_form} (which is approximately {numeric_value})")

solve_equation()
