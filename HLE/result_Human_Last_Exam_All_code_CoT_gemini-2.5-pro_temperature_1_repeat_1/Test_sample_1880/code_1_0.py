import sympy
from sympy import symbols, Poly, solve, sqrt

def solve_galois_order():
    """
    Computes the order of the Galois group for the polynomial x^4 + 8x + 14.
    """
    x, y = symbols('x y')

    # Step 1: Define the polynomial and check irreducibility
    p = 8
    q = 14
    f = x**4 + p*x + q
    print(f"The polynomial is f(x) = x^4 + {p}x + {q}")

    # For p=2, the prime 2 divides 14 and 8, but not 1. And 2^2=4 does not divide 14.
    # By Eisenstein's criterion with p=2, the polynomial is irreducible over Q.
    # We can also verify with sympy.
    is_irred = Poly(f, x, domain='QQ').is_irreducible
    print(f"Is the polynomial irreducible over Q? {is_irred}\n")
    if not is_irred:
        print("The polynomial is reducible, the method for irreducible quartics doesn't apply directly.")
        return

    # Step 2: Form the resolvent cubic
    # For x^4 + px + q, the resolvent is y^3 - 4qy - p^2
    g_p = -4 * q
    g_q = -p**2
    g = y**3 + g_p*y + g_q
    print(f"The resolvent cubic is g(y) = y^3 + ({g_p})y + ({g_q})")

    # Step 3: Find and analyze the roots of the resolvent cubic
    resolvent_roots = solve(g, y)
    print(f"The roots of the resolvent cubic are: {resolvent_roots}\n")

    rational_roots = [r for r in resolvent_roots if r.is_rational]
    num_rational_roots = len(rational_roots)
    print(f"Number of rational roots of the resolvent: {num_rational_roots}")

    # Step 4: Distinguish between candidate groups
    # Since there is exactly one rational root, the Galois group is either D_4 or C_4.
    # The group is C_4 if f(x) is reducible over the field K defined by the
    # other roots of the resolvent, and D_4 otherwise.
    # The field K is Q(sqrt(D)) where D is the discriminant of the quadratic factor of g(y).
    # The irrational roots are -4 +/- 2*sqrt(2), so the field is Q(sqrt(2)).
    # The test for reducibility is to check if any root of the resolvent is a square in K.
    
    print("The resolvent has one rational root, so the Galois group is a subgroup of D_4.")
    print("It is either D_4 (order 8) or C_4 (order 4).")
    
    # The quadratic field is Q(sqrt(2))
    K_sqrt_term = sqrt(2)
    print(f"The quadratic field K is Q({K_sqrt_term})\n")
    
    print("To distinguish between D_4 and C_4, we check if any root of the resolvent cubic is a square in K.")
    
    is_reducible_over_K = False
    for r in resolvent_roots:
        # Let's check if sqrt(r) can be written as a + b*sqrt(2) for rational a, b.
        # Let sqrt(r) = a + b*K_sqrt_term. Then r = (a^2 + 2b^2) + 2ab*K_sqrt_term.
        # We solve for rational a, b.

        # Case 1: r is rational. r = r_rat + 0*K_sqrt_term.
        if r.is_rational:
            # We need 2ab = 0 and a^2 + 2b^2 = r.
            # If a=0: 2b^2 = r => b = sqrt(r/2). Rational if r/2 is a square of a rational.
            # If b=0: a^2 = r => a = sqrt(r). Rational if r is a square of a rational.
            if sqrt(r).is_rational or sqrt(r/2).is_rational:
                is_reducible_over_K = True
                print(f"Checking root r = {r}:")
                print(f"  sqrt({r}) = {sympy.simplify(sqrt(r))}, which is in the field K = Q({K_sqrt_term}).")
                break
        
        # Case 2: r is irrational, of the form c + d*K_sqrt_term
        else:
            # r = c + d*K_sqrt_term = (a^2 + 2b^2) + 2ab*K_sqrt_term
            c, d = r.as_coeff_add(K_sqrt_term)
            d = d[0][0] # Extract coefficient
            # We need 2ab = d and a^2 + 2b^2 = c.
            # b = d/(2a) => a^2 + 2*(d/(2a))^2 = c => a^2 + d^2/(2a^2) = c
            # 2a^4 - 2ca^2 + d^2 = 0. Let u = a^2. 2u^2 - 2cu + d^2 = 0
            u_coeffs = [2, -2*c, d**2]
            u_sols = solve(u_coeffs[0]*x**2 + u_coeffs[1]*x + u_coeffs[2], x)
            # We need a rational solution for u=a^2, where u is a positive rational whose sqrt is rational.
            found_sol = False
            for u in u_sols:
                if u.is_rational and u >= 0 and sqrt(u).is_rational:
                    found_sol = True
                    break
            if found_sol:
                is_reducible_over_K = True
                print(f"Checking root r = {r}: It is a square in K.")
                break
            else:
                 print(f"Checking root r = {r}: It is not a square in K.")

    # Step 5: Final conclusion
    print("\nConclusion:")
    if is_reducible_over_K:
        print("The polynomial is reducible over K, so the Galois group is C_4 (cyclic group of order 4).")
        order = 4
    else:
        print("The polynomial is irreducible over K, so the Galois group is D_4 (dihedral group of order 8).")
        order = 8

    print(f"\nThe order of the Galois group is {order}.")
    return order

if __name__ == '__main__':
    solve_galois_order()