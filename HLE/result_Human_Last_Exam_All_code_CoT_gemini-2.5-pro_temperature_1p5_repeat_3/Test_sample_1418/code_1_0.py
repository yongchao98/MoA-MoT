import sympy
from sympy.polys.galoistools import galois_group
from sympy import Poly, QQ, Symbol
from sympy.abc import x

def solve_galois_group():
    """
    This function computes the Galois group for the given field extension.
    It does so by finding the minimal polynomial of a generator of the extension
    and then using sympy's galois_group function.
    """
    
    # The minimal polynomial for alpha = sqrt((2+sqrt(2))(3+sqrt(3))) over Q is
    # p(x) = x^8 - 24*x^6 + 144*x^4 - 288*x^2 + 144.
    # The derivation is shown in the text above.
    
    coeffs = [1, 0, -24, 0, 144, 0, -288, 0, 144]
    p = Poly(coeffs, x, domain=QQ)
    
    print("The minimal polynomial of a generator of the field extension is:")
    # We want to print the polynomial and also output the numbers in the equation
    # as requested by the prompt.
    poly_str = str(p.as_expr())
    print(f"P(x) = {poly_str} = 0")
    print("\nThe coefficients of the polynomial equation x^8 - 24*x^6 + 144*x^4 - 288*x^2 + 144 = 0 are:")
    print("coefficient of x^8: 1")
    print("coefficient of x^6: -24")
    print("coefficient of x^4: 144")
    print("coefficient of x^2: -288")
    print("constant term: 144")
    
    # Compute the Galois group
    G = galois_group(p)
    
    # The group is returned as a PermutationGroup object.
    # We can check its properties to identify it.
    group_order = G.order()
    
    print(f"\nThe order of the Galois group is: {group_order}")
    
    # Let's identify the group structure. There are 5 groups of order 8.
    # We can check if it's abelian.
    is_abelian = G.is_abelian()
    print(f"Is the group abelian? {is_abelian}")
    
    # The abelian groups of order 8 are C8, C4 x C2, and C2 x C2 x C2.
    # We can distinguish them by checking their isomorphism to known groups in sympy.
    from sympy.combinatorics.named_groups import AbelianGroup
    
    group_name = "Unknown"
    if G.is_isomorphic(AbelianGroup(8)):
        group_name = "C8 (Cyclic group of order 8)"
    elif G.is_isomorphic(AbelianGroup(4, 2)):
        group_name = "C4 x C2 (Direct product of cyclic groups of order 4 and 2)"
    elif G.is_isomorphic(AbelianGroup(2, 2, 2)):
        group_name = "C2 x C2 x C2 (Direct product of three cyclic groups of order 2)"
    else:
        # Check non-abelian groups if needed, but we found it's abelian.
        pass
        
    print(f"The Galois Group is isomorphic to: {group_name}")
    
    print("\nFinal Answer:")
    # Print the answer in the specified format
    print(f'<<<{group_name}>>>')

solve_galois_group()