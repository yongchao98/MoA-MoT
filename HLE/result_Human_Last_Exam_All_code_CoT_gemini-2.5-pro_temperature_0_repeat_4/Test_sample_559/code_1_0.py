import sympy as sp

def find_separatrix_equation():
    """
    This function finds the equation of the separatrix for the given system of ODEs.

    The system is:
    d'(t) = 2*d**2 + (5*u**2 - 3*u)*d - u**3*(1-u)
    u'(t) = u**2 * (u - 1)

    The method involves finding a particular solution to the trajectory equation
    dd/du = d_prime / u_prime of the form d = c*u**2.
    """
    # Define the symbolic variables for u and the constant c
    u, c = sp.symbols('u c')

    # 1. Define the hypothesized form of the separatrix and its derivative w.r.t. u
    # We guess the separatrix is of the form d = c * u^2
    d_guess = c * u**2
    dd_du_guess = sp.diff(d_guess, u)
    k_val = 2

    # 2. Define the original ODEs, d'(t) and u'(t), in terms of u and our guess for d
    d_prime = 2*d_guess**2 + (5*u**2 - 3*u)*d_guess - u**3*(1-u)
    u_prime = u**2 * (u - 1)

    # 3. Form the trajectory equation: dd/du = d_prime / u_prime
    # We create an algebraic equation that must be true if our guess is correct.
    # dd_du_guess * u_prime - d_prime = 0
    equation = sp.simplify(dd_du_guess * u_prime - d_prime)

    # 4. For the equation to hold for all u, it must be a zero polynomial in u.
    # We extract the coefficients of the polynomial in u and set them to zero.
    poly_in_u = sp.Poly(equation, u)
    coefficients = poly_in_u.coeffs()

    # 5. Solve for the constant c using the system of coefficient equations.
    # The `solve` function can find the value of c that makes all coefficients zero.
    solutions_for_c = sp.solve(coefficients, c)

    # 6. Print the results.
    if solutions_for_c:
        # We expect a single constant solution for c
        c_val = solutions_for_c[0]
        
        print("The separatrix is an invariant curve, which is a solution to the system's trajectory equation.")
        print(f"We searched for a solution of the form: d = c * u^k")
        print("\nBy substituting this form into the equation and solving, we find the constants:")
        print(f"The power k is: {k_val}")
        print(f"The coefficient c is: {c_val}")
        
        # Construct the final equation string
        final_equation = c_val * u**k_val
        print("\nThus, the equation for the separatrix is:")
        print(f"d = {sp.pretty(final_equation, use_unicode=True)}")

    else:
        print("A solution of the form d = c*u^2 could not be found.")

# Run the function to find and print the separatrix equation
find_separatrix_equation()