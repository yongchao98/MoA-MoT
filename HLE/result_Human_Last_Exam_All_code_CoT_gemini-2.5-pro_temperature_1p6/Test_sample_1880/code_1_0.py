import sympy
from sympy import symbols, Poly, Rational, sqrt, factor

def compute_galois_group_order():
    """
    Computes the order of the Galois group for the polynomial x^4 + 8x + 14.
    """
    # Step 1: Define the polynomial and its parameters for the depressed quartic form x^4 + px + q
    x = symbols('x')
    f_poly = x**4 + 8*x + 14
    p = 8
    q = 14
    
    print(f"Analyzing the polynomial f(x) = {f_poly}")
    
    # Check for irreducibility using sympy (or manually by Eisenstein's criterion with p=2)
    is_irred = Poly(f_poly, x, domain='QQ').is_irreducible()
    if not is_irred:
        print("The polynomial is reducible, the analysis for irreducible quartics doesn't apply directly.")
        return

    print("Step 1: The polynomial is irreducible over the rational numbers Q.")
    
    # Step 2: Calculate the discriminant of the quartic
    # For x^4 + px + q, the discriminant is 256*q^3 - 27*p^4
    discriminant = 256 * q**3 - 27 * p**4
    print(f"\nStep 2: The discriminant is calculated as 256*({q})^3 - 27*({p})^4 = {discriminant}")
    
    # Check if the discriminant is a perfect square
    sqrt_discriminant = sympy.integer_nthroot(discriminant, 2)
    is_square = sqrt_discriminant[1]
    
    if is_square:
        print(f"The discriminant {discriminant} is a perfect square. The Galois group is a subgroup of A_4.")
    else:
        print(f"The discriminant {discriminant} is not a perfect square. The Galois group is not a subgroup of A_4.")
        print("Possible groups are S_4, D_4, or C_4.")
        
    # Step 3: Form the resolvent cubic
    # For x^4 + px + q, the resolvent cubic is y^3 - 4*q*y - p^2
    y = symbols('y')
    resolvent_cubic = y**3 - 4*q*y - p**2
    print(f"\nStep 3: The resolvent cubic is g(y) = y^3 - 4*({q})*y - ({p})^2 = {resolvent_cubic}")
    
    # Step 4: Analyze the roots of the resolvent cubic
    resolvent_roots = sympy.roots(resolvent_cubic, y)
    rational_roots = [r for r in resolvent_roots if r.is_rational]
    num_rational_roots = len(rational_roots)
    
    print(f"\nStep 4: The resolvent cubic has {num_rational_roots} rational root(s).")
    
    # Step 5: Identify the Galois Group
    galois_group_name = ""
    order = 0
    
    print("\nStep 5: Identifying the Galois group...")
    if num_rational_roots == 0:
        # Resolvent is irreducible
        if is_square:
             # This case is A_4
             galois_group_name = "A_4"
             order = 12
        else:
             # This case is S_4
             galois_group_name = "S_4"
             order = 24
        print(f"The resolvent is irreducible. Combined with the discriminant info, the group is {galois_group_name}.")

    elif num_rational_roots >= 1: # Reducible case
        if num_rational_roots == 3:
            # This case is V_4
            galois_group_name = "V_4"
            order = 4
            print(f"The resolvent has 3 rational roots. The group is {galois_group_name}.")
        else: # num_rational_roots == 1
            print("The resolvent is reducible with one rational root. The group is either D_4 or C_4.")
            # Distinguish D_4 and C_4
            r1 = rational_roots[0]
            quad_factor = sympy.div(resolvent_cubic, y - r1)[0]
            # Discriminant of the quadratic factor
            a, b, c = Poly(quad_factor, y).all_coeffs()
            quad_discriminant = b**2 - 4*a*c
            
            # Check if f(x) is reducible over Q(sqrt(quad_discriminant))
            is_reducible_over_ext = not Poly(f_poly, x, extension=sqrt(quad_discriminant)).is_irreducible()
        
            if is_reducible_over_ext:
                print(f"The polynomial f(x) is reducible over the field extension defined by the other resolvent roots.")
                galois_group_name = "D_4"
                order = 8
            else:
                print(f"The polynomial f(x) is irreducible over the field extension defined by the other resolvent roots.")
                galois_group_name = "C_4"
                order = 4
            print(f"The final test shows the Galois group is {galois_group_name}.")

    print(f"\nConclusion: The Galois group for {f_poly} is {galois_group_name}.")
    print(f"The order of this group is {order}.")


compute_galois_group_order()