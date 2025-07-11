import numpy as np
from sympy import Poly, symbols

def get_invariant_poly(dim, name):
    """
    Computes the Poincare polynomial of invariants for a derivation
    that acts as a single nilpotent Jordan block on a vector space of dimension `dim`.
    The dimension of invariants at each degree q is known to be 1.
    """
    x = symbols('x')
    coeffs = [1] * (dim + 1)
    poly = Poly(list(reversed(coeffs)), x)
    print(f"Subsystem {name}:")
    print(f" The space is {dim}-dimensional.")
    print(f" The Poincare polynomial of its invariants is: {poly.as_expr()}")
    print("-" * 20)
    return poly

def main():
    """
    Computes the Poincare polynomial of the given 6-dimensional Lie algebra.
    """
    print("Step 1: Decompose the problem.")
    print("The Lie algebra g is a semidirect product R x h, where h is a 5D abelian ideal.")
    print("Its Poincare polynomial is P_g(x) = (1+x) * P_inv(x).")
    print("The action on the dual space h* splits into two blocks: a 3D block and a 2D block.")
    print("Thus, P_inv(x) = P_inv,1(x) * P_inv,2(x).")
    print("-" * 20)

    # Step 2: Compute Poincare polynomial for the first invariant subsystem
    # This corresponds to the subspace spanned by {e^2, e^3, e^4}
    # with derivation d(e^2)=-e^3, d(e^3)=-e^4, d(e^4)=0.
    # This is a nilpotent block of size 3.
    p_inv1 = get_invariant_poly(dim=3, name="1")

    # Step 3: Compute Poincare polynomial for the second invariant subsystem
    # This corresponds to the subspace spanned by {e^5, e^6}
    # with derivation d(e^5)=-e^6, d(e^6)=0.
    # This is a nilpotent block of size 2.
    p_inv2 = get_invariant_poly(dim=2, name="2")

    # Step 4: Combine the polynomials
    print("Step 2: Combine the results.")
    p_inv = p_inv1 * p_inv2
    print(f"The total Poincare polynomial of invariants is:")
    print(f"P_inv(x) = P_inv,1(x) * P_inv,2(x) = {p_inv.as_expr()}")
    print("-" * 20)

    x = symbols('x')
    one_plus_x = Poly(1 + x, x)
    p_g = p_inv * one_plus_x

    # Step 5: Print the final result
    print("Step 3: Compute the final Poincare polynomial for g.")
    print(f"P_g(x) = (1+x) * P_inv(x)")
    
    # Format the final polynomial string
    final_poly_str = str(p_g.as_expr()).replace('**', '^')
    
    print("\nThe final Poincare polynomial is:")
    print(final_poly_str)

    # Also printing the equation with each coefficient as requested
    coeffs = p_g.all_coeffs()
    terms = []
    degree = p_g.degree()
    for i, c in enumerate(coeffs):
        if c == 0:
            continue
        power = degree - i
        if power == 0:
            terms.append(str(c))
        elif power == 1:
            terms.append(f"{c}*x")
        else:
            terms.append(f"{c}*x^{power}")
    
    equation = " + ".join(terms).replace("+ -", "- ")
    
    print("\nThe polynomial written as an equation is:")
    print(equation)


if __name__ == '__main__':
    main()
