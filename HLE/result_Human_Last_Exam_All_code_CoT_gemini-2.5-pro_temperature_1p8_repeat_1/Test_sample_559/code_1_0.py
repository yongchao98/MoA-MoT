import sympy as sp

def find_separatrix():
    """
    This function finds the separatrix for the given system of differential equations
    by systematically identifying the saddle point and its unstable manifold.
    """
    # Define the variables and the system of differential equations
    u, d = sp.symbols('u d')
    u_prime = (u - 1) * u**2
    d_prime = 2*d**2 + (-3*u + 5*u**2)*d - u*(1 - u)*u**2

    # Find equilibrium points by solving d'(t)=0 and u'(t)=0
    equilibria = sp.solve([u_prime, d_prime], [u, d])
    
    # Compute the Jacobian matrix of the system
    J = sp.Matrix([u_prime, d_prime]).jacobian([u, d])

    # Identify the saddle point (1, -1)
    saddle_point = (1, -1)
    J_saddle = J.subs({u: saddle_point[0], d: saddle_point[1]})

    # Find the eigenvectors of the Jacobian at the saddle point
    eigenvects = J_saddle.eigenvects()
    
    unstable_eigenvector = None
    # The unstable manifold corresponds to the positive eigenvalue
    for eigval, multiplicity, vectors in eigenvects:
        if eigval > 0:
            unstable_eigenvector = vectors[0]
            break

    # Calculate the slope of the unstable manifold at the saddle point
    m_saddle = unstable_eigenvector[1] / unstable_eigenvector[0]

    # Assume the separatrix is a quadratic polynomial d = a*u**2 + b*u + c
    a, b, c = sp.symbols('a b c')
    d_poly = a*u**2 + b*u + c

    # Set up equations to find the coefficients a, b, c
    # The separatrix connects the saddle (1,-1) to the equilibrium (0,0)
    # Condition 1: The curve passes through (0,0)
    eq1 = d_poly.subs(u, 0) - 0
    # Condition 2: The curve passes through the saddle point (1,-1)
    eq2 = d_poly.subs(u, saddle_point[0]) - saddle_point[1]
    # Condition 3: The slope at the saddle point matches the unstable eigenvector
    d_poly_derivative = sp.diff(d_poly, u)
    eq3 = d_poly_derivative.subs(u, saddle_point[0]) - m_saddle
    
    # Solve for the coefficients
    coeffs = sp.solve([eq1, eq2, eq3], [a, b, c])
    a_val = coeffs[a]
    b_val = coeffs[b]
    c_val = coeffs[c]

    # Print the equation of the separatrix and its constituent numbers
    print("The equation of the separatrix is d = a * u**2 + b * u + c.")
    print("The numbers in this equation are:")
    print(f"a = {a_val}")
    print(f"b = {b_val}")
    print(f"c = {c_val}")
    print("\nThe simplified equation is d = -u**2, or d + u**2 = 0.")

find_separatrix()