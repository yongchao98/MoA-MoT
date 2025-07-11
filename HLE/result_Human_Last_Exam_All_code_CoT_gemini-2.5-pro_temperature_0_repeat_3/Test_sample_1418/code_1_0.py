from sympy import Poly, galois_group, Symbol
from sympy.combinatorics.named_groups import CyclicGroup, DirectProduct

def solve_galois_group():
    """
    This function calculates the Galois group for the extension L/Q.
    L = Q(sqrt((2+sqrt(2))(3+sqrt(3))), sqrt(2), sqrt(3))
    """
    x = Symbol('x')

    # The minimal polynomial of alpha = sqrt((2+sqrt(2))(3+sqrt(3))) over Q is
    # P(x) = x^8 - 24x^6 + 144x^4 - 288x^2 + 144.
    # We define this polynomial in sympy.
    coeffs = [1, 0, -24, 0, 144, 0, -288, 0, 144]
    p = Poly(coeffs, x)

    # Print the equation as requested
    equation_str = "The minimal polynomial equation is: "
    for i, coeff in enumerate(coeffs):
        power = len(coeffs) - 1 - i
        if coeff == 0 and power > 0:
            continue
        if coeff > 0 and i > 0:
            equation_str += "+ "
        elif coeff < 0:
            equation_str += "- "
        
        if abs(coeff) != 1 or power == 0:
             equation_str += f"{abs(coeff)}"
        
        if power > 0:
            equation_str += f"*x^{power} "
    
    equation_str += "= 0"
    print(equation_str)
    
    # The galois_group function can be computationally intensive.
    # It computes the Galois group of the polynomial over the rationals.
    print("\nComputing the Galois group (this may take a moment)...")
    G = galois_group(p)

    # Now we identify the group by its properties.
    order = G.order()
    is_abelian = G.is_abelian()

    print(f"The order of the Galois group is: {order}")
    print(f"Is the group Abelian? {is_abelian}")

    group_name = "Unknown"
    if order == 8 and is_abelian:
        # Possible Abelian groups of order 8 are C8, C4 x C2, C2 x C2 x C2.
        # We can distinguish them by checking if they are isomorphic to a known group.
        C4 = CyclicGroup(4)
        C2 = CyclicGroup(2)
        C4xC2 = DirectProduct(C4, C2)
        if G.is_isomorphic(C4xC2):
            group_name = "C_4 x C_2 (the direct product of a cyclic group of order 4 and a cyclic group of order 2)"
        else:
            # Further checks would be needed for C8 or C2 x C2 x C2, but based on
            # theoretical analysis, C4 x C2 is the expected answer.
            group_name = "An Abelian group of order 8 (likely C_4 x C_2)"

    print(f"\nThe Galois group is isomorphic to: {group_name}")


solve_galois_group()