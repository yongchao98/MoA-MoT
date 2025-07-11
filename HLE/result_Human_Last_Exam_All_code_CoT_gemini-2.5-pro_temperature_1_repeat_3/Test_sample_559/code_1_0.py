import sympy as sp

def find_separatrix():
    """
    This function finds the separatrix for the given system of differential equations
    by assuming a polynomial form and solving for its coefficients.
    """
    # Step 1: Define symbols and the proposed form of the separatrix.
    # We are working in the (u, d) phase plane.
    # d is treated as a function of u along the trajectory.
    u, a, d = sp.symbols('u a d')

    # We hypothesize that the separatrix is a curve of the form d = a * u^k.
    # Balancing the powers in the ODEs suggests k=2 is a good candidate.
    # So we look for a separatrix of the form d = a * u**2.
    d_guess = a * u**2
    print("Searching for a separatrix of the form: d = a * u**2")
    
    # Step 2: Define the system of differential equations from the problem statement.
    # d'(t) from the first ODE of the system.
    d_dot_system = 2 * d**2 + (-3*u + 5*u**2)*d - u**3*(1-u)
    
    # u'(t) from the second ODE of the system.
    u_dot_system = u**2 * (u - 1)

    # Step 3: Apply the condition for an invariant manifold.
    # For a curve d(u) to be a trajectory, the flow vector (u', d') must be tangent to it.
    # This means d'(t) must equal (d(d)/du) * u'(t).
    
    # First, calculate d'(t) by substituting our guess d = a*u**2 into the first ODE.
    d_dot_on_curve = d_dot_system.subs(d, d_guess)

    # Second, calculate d'(t) using the chain rule on our guess d = a*u**2,
    # and substituting the expression for u'(t) from the second ODE.
    d_dot_from_chain_rule = sp.diff(d_guess, u) * u_dot_system

    # Step 4: The two expressions for d'(t) must be equal.
    # Their difference must be zero for all values of u.
    equation = sp.simplify(d_dot_on_curve - d_dot_from_chain_rule)

    # Step 5: The resulting equation is a polynomial in u. For it to be zero for all u,
    # all its coefficients must be zero. We can extract these coefficients.
    # The `as_poly()` method is convenient for this. If the expression is not a polynomial
    # or is zero, it returns None.
    if equation == 0:
        print("The guess is correct for any 'a', which is unlikely. Check system.")
        return
        
    poly_equation = sp.Poly(equation, u)
    coefficients = poly_equation.coeffs()

    # Step 6: Create algebraic equations for 'a' by setting the coefficients to zero.
    # We can solve this system of equations for 'a'. Usually, the simplest equation is enough.
    # Let's use the last non-zero coefficient, which corresponds to the lowest power of u.
    # From manual derivation: equation = u**4 * (-2*a**2 - 3*a - 1) + u**3 * (a + 1)
    # The coefficient of u^3 gives a + 1 = 0, which is simple to solve.
    coeff_eq_for_a = coefficients[-1]
    solutions_for_a = sp.solve(coeff_eq_for_a, a)
    
    # We expect a single solution for 'a'.
    if not solutions_for_a:
        print("Could not find a value for 'a'. The separatrix may have a different form.")
        return
        
    a_val = solutions_for_a[0]
    
    # We can verify this solution in the other coefficient equations.
    # For example, -2*a**2 - 3*a - 1 should be 0 for a = -1.
    # -2*(-1)**2 - 3*(-1) - 1 = -2 + 3 - 1 = 0. The value is consistent.

    print(f"\nSolved for the coefficient 'a', found a = {a_val}")

    # Step 7: Construct and print the final equation of the separatrix.
    # The equation is d = a_val * u**2, which we can write as 1*d - a_val*u**2 = 0.
    c_d = 1
    c_u2 = -a_val

    print("\nThe separatrix is given by the equation d = -u^2.")
    print("To meet the output format requirement, we write it as:")
    print(f"{c_d}*d + {c_u2}*u**2 = 0")


if __name__ == '__main__':
    find_separatrix()