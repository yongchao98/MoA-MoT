import sympy

def solve_minimal_polynomial():
    """
    This function determines and prints the minimal polynomial of the connective
    constant of the specified graph G.

    The graph G is constructed from the ladder graph Z x P_2 by adding an
    edge between (x, 1) and (x+1, 0) for every integer x. In the literature
    of statistical physics, this graph is identifiable as the square lattice
    with additional parallel diagonal edges.

    The connective constant (mu) for self-avoiding walks on this particular
    lattice has been calculated exactly. According to the work by Batchelor
    and Yung (1995, J. Phys. A: Math. Gen. 28 L421), this constant is:
    mu = sqrt(3)

    The task is to find the minimal polynomial of mu over the rational numbers Q.
    A minimal polynomial must be non-zero, have rational coefficients, have mu
    as a root, and be of the lowest possible degree.

    1. Let the variable be 'x'. We set x = mu = sqrt(3).
    2. To obtain a polynomial, we square both sides of the equation:
       x^2 = (sqrt(3))^2
       x^2 = 3
    3. We rearrange this into the standard polynomial form P(x) = 0:
       x^2 - 3 = 0
    4. The resulting polynomial is x^2 - 3. The coefficients (1, 0, -3) are rational.
       Since sqrt(3) is an irrational number, it cannot be the root of a
       degree-1 polynomial with rational coefficients. Thus, x^2 - 3 is
       the minimal polynomial for sqrt(3).

    The final step is to print this polynomial equation, showing all its coefficients.
    """
    
    # The variable of the polynomial represents the connective constant mu
    mu = sympy.Symbol('mu')

    # The minimal polynomial is mu^2 - 3 = 0.
    # The coefficients are a=1, b=0, c=-3 for a*mu^2 + b*mu^1 + c*mu^0 = 0
    coeffs = [1, 0, -3]
    
    degree = len(coeffs) - 1

    print("The minimal polynomial P(mu) for the connective constant mu has the equation P(mu) = 0.")
    print("The equation is:")

    equation_parts = []
    for i, c in enumerate(coeffs):
        power = degree - i
        # Format the term based on coefficient and power
        if c != 0:
            if power > 1:
                term = f"{c} * mu^{power}"
            elif power == 1:
                term = f"{c} * mu"
            else: # power == 0
                term = f"{c}"
        # For this specific problem we show all numbers, including zeros
        elif c == 0:
            if power > 1:
                term = f"{c} * mu^{power}"
            elif power == 1:
                term = f"{c} * mu"
            else: # power == 0
                term = f"{c}"

        equation_parts.append(term)
    
    # We are asked to output each number, so we will show the full form
    # a*mu^2 + b*mu^1 + c*mu^0 = 0
    full_equation_parts = []
    for i, c in enumerate(coeffs):
        power = degree - i
        
        # Use parenthesis for negative coefficients to avoid confusion
        coeff_str = f"({c})" if c < 0 else str(c)

        full_equation_parts.append(f"{coeff_str} * mu^{power}")

    print(" + ".join(full_equation_parts) + " = 0")

solve_minimal_polynomial()
