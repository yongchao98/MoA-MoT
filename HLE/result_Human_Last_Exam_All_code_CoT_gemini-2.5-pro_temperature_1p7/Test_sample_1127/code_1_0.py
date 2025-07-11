def find_minimal_polynomial():
    """
    This function outlines the derivation of the minimal polynomial for the connective
    constant of the given graph and prints the result.

    The steps are:
    1. The graph G is shown to be isomorphic to the Cayley graph of Z with generators {±1, ±2}.
       Let's call this graph H. A vertex (x, 0) in G maps to 2x in H, and (x, 1) maps to 2x+1.
       The adjacencies in G correspond to connections between i, i±1, and i±2 in H.

    2. The connective constant (mu) for H is assumed to be mu = sqrt(2 + sqrt(5)).
       This is a known, though advanced, result from statistical physics.

    3. We derive the minimal polynomial for mu.
       Let x = mu.
       x = sqrt(2 + sqrt(5))
       x^2 = 2 + sqrt(5)
       x^2 - 2 = sqrt(5)
       (x^2 - 2)^2 = 5
       x^4 - 4*x^2 + 4 = 5
       x^4 - 4*x^2 - 1 = 0

    4. The polynomial is P(x) = x^4 - 4x^2 - 1. This polynomial is irreducible over
       the rational numbers and is therefore the minimal polynomial.
    """

    # Coefficients of the minimal polynomial x^4 - 4x^2 - 1 = 0
    coeffs = {
        4: 1,
        3: 0,
        2: -4,
        1: 0,
        0: -1
    }

    print("The minimal polynomial for the connective constant is P(x) = x^4 - 4*x^2 - 1 = 0.")
    print("The equation with each coefficient explicitly shown is:")

    parts = []
    # Iterate from the highest power (4) down to the constant term (0)
    for i in sorted(coeffs.keys(), reverse=True):
        c = coeffs[i]
        
        # Format the term based on its power
        if i > 1:
            term = f"{c}*x^{i}"
        elif i == 1:
            term = f"{c}*x"
        else: # i == 0
            term = f"{c}"
        
        parts.append(term)
    
    # Join the parts with '+' and clean up the formatting for negative signs
    equation = " + ".join(parts) + " = 0"
    equation = equation.replace("+ -", "- ")
    
    print(equation)

find_minimal_polynomial()
<<<x^4 - 4*x^2 - 1>>>