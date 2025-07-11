import sympy

def find_and_analyze_separatrix():
    """
    Finds and verifies the separatrix for the given system of ODEs.
    """
    # Define symbols for the variables
    d, u = sympy.symbols('d u')

    # Define the system of differential equations from the problem
    d_prime = 2*d**2 + (-3*u + 5*u**2)*d - u**3*(1 - u)
    u_prime = (u - 1)*u**2

    print("System of differential equations:")
    print(f"d'(t) = {d_prime}")
    print(f"u'(t) = {u_prime}\n")

    # Step 1: Find equilibrium points by solving d'=0 and u'=0
    print("Step 1: Finding equilibrium points...")
    try:
        equilibria = sympy.solve([d_prime, u_prime], (d, u), dict=True)
        print(f"The equilibrium points (d, u) are: {[(p[d], p[u]) for p in equilibria]}\n")
    except Exception as e:
        print(f"Could not solve for equilibrium points symbolically: {e}")
        return

    # Step 2: Classify equilibrium points using the Jacobian matrix
    print("Step 2: Classifying equilibrium points...")
    J = sympy.Matrix([d_prime, u_prime]).jacobian([d, u])
    
    saddle_points = []
    for point in equilibria:
        d_val, u_val = point[d], point[u]
        # We are interested in the u >= 0 half-plane
        if u_val < 0:
            continue
            
        print(f"--- Analyzing point ({d_val}, {u_val}) ---")
        J_at_point = J.subs(point)
        
        # Check for degenerate cases where linearization fails
        if J_at_point.is_zero_matrix:
            print("Jacobian is the zero matrix. Linearization fails. This is a non-hyperbolic point.")
            continue
            
        eigenvals = J_at_point.eigenvals()
        print(f"Jacobian eigenvalues: {list(eigenvals.keys())}")
        
        # Classify based on eigenvalues
        real_parts = [sympy.re(ev) for ev in eigenvals]
        if any(r > 0 for r in real_parts) and any(r < 0 for r in real_parts):
            print("Classification: Saddle point")
            saddle_points.append(point)
        elif all(r > 0 for r in real_parts):
            print("Classification: Unstable node (source)")
        elif all(r < 0 for r in real_parts):
            print("Classification: Stable node (sink)")
        else:
            print("Classification: Non-hyperbolic (center, spiral, or other)")
    print("---\n")

    if not saddle_points:
        print("No saddle points found. Cannot determine separatrix using this method.")
        return

    # Step 3: Propose and verify the separatrix equation
    print("Step 3: Verifying the proposed separatrix d = -u^2...")
    
    # The proposed separatrix is d = -u^2
    separatrix_d = -u**2
    
    # An invariant curve must satisfy d'(t) = (dd/du) * u'(t)
    # Let's verify this identity.
    
    # Left-hand side: substitute d = -u^2 into d'
    lhs = d_prime.subs(d, separatrix_d)
    lhs_simplified = sympy.simplify(lhs)
    print(f"Substituting d = -u^2 into d'(t) gives: {lhs_simplified}")

    # Right-hand side: (dd/du) * u'
    d_separatrix_du = sympy.diff(separatrix_d, u)
    rhs = d_separatrix_du * u_prime
    rhs_simplified = sympy.simplify(rhs)
    print(f"Calculating (dd/du) * u'(t) gives: {rhs_simplified}")
    
    # Check if LHS equals RHS
    if sympy.simplify(lhs - rhs) == 0:
        print("Verification successful: The curve d = -u^2 is an invariant curve.\n")
    else:
        print("Verification failed: The curve d = -u^2 is not an invariant curve.\n")
        return

    # Step 4: Check if the separatrix passes through a saddle point
    print("Step 4: Checking if the invariant curve passes through a saddle point...")
    saddle_point = saddle_points[0] # Using the first found saddle
    d_saddle, u_saddle = saddle_point[d], saddle_point[u]
    
    check_val = separatrix_d.subs(u, u_saddle)
    if check_val == d_saddle:
        print(f"The invariant curve d = -u^2 passes through the saddle point ({d_saddle}, {u_saddle}).")
        print("Therefore, it is a separatrix of the system.\n")
    else:
        print(f"The invariant curve does not pass through the saddle point ({d_saddle}, {u_saddle}).")
        return

    # Final result
    print("The final equation for the separatrix is:")
    # We want to print each number in the equation, so we make coefficients explicit.
    coeff_d = 1
    coeff_u2 = 1
    # Representing as d + u^2 = 0
    print(f"{coeff_d}*d + {coeff_u2}*u**2 = 0")
    # Or as d = -u^2
    coeff_u2_rhs = -1
    print(f"Alternatively: d = {coeff_u2_rhs} * u**2")


if __name__ == "__main__":
    find_and_analyze_separatrix()