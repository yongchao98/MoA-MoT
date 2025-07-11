import sympy

def find_separatrix():
    """
    This function analyzes a system of ODEs to find its separatrix.
    It follows these steps:
    1. Defines the system of ODEs symbolically.
    2. Solves for the equilibrium points.
    3. Performs a stability analysis for each equilibrium point using the Jacobian matrix.
    4. Identifies the saddle point.
    5. Tests a candidate curve and verifies it is an invariant manifold passing
       through the saddle point, which makes it a separatrix.
    6. Prints the final equation for the separatrix.
    """
    # Step 1: Define variables and the system of ODEs
    d, u = sympy.symbols('d u')
    d_prime = 2*d**2 + (-3*u + 5*u**2)*d - u**3*(1-u)
    u_prime = u**2 * (u-1)

    print("--- System of Differential Equations ---")
    print(f"d'(t) = {d_prime}")
    print(f"u'(t) = {u_prime}\n")

    # Step 2: Find equilibrium points
    # Solve d'(t)=0 and u'(t)=0 for (d, u)
    try:
        equilibria = sympy.solve([d_prime, u_prime], (d, u), dict=True)
    except NotImplementedError: # Handle cases where sympy.solve might fail
        print("Could not solve for equilibria automatically. Analyzing manually.")
        equilibria = [{d: 0, u: 0}, {d: 0, u: 1}, {d: -1, u: 1}]

    print("--- Equilibrium Points (d, u) ---")
    for point in equilibria:
        print(f"({point[d]}, {point[u]})")
    print("")

    # Step 3: Linearization and Stability Analysis
    J = sympy.Matrix([d_prime, u_prime]).jacobian([d, u])
    
    print("--- Stability Analysis ---")
    saddle_point = None
    for point in equilibria:
        d_val, u_val = point[d], point[u]
        J_at_point = J.subs(point)
        eigenvals = J_at_point.eigenvals()
        
        real_parts = [sympy.re(ev) for ev in eigenvals.keys()]
        
        stability = "Undetermined"
        if J_at_point == sympy.zeros(2):
            stability = "Non-hyperbolic (Linearization is inconclusive)"
        elif any(r > 0 for r in real_parts) and any(r < 0 for r in real_parts):
            stability = "Saddle Point"
            saddle_point = point
        elif all(r > 0 for r in real_parts):
            stability = "Unstable Node"
        elif all(r < 0 for r in real_parts):
            stability = "Stable Node"

        print(f"Point ({d_val}, {u_val}):")
        print(f"  Eigenvalues: {list(eigenvals.keys())}")
        print(f"  Stability Type: {stability}\n")

    # Step 4 & 5: Find and verify the separatrix equation
    print("--- Finding the Separatrix ---")
    if not saddle_point:
        print("No saddle point found. Cannot determine the primary separatrix in this manner.")
        return

    print(f"The saddle point is at ({saddle_point[d]}, {saddle_point[u]}).")
    print("Separatrices are often the stable/unstable manifolds of saddle points.")
    print("We hypothesize a solution of the form d = -u^2 and verify if it's a trajectory.")

    # Candidate separatrix equation
    separatrix_candidate = -u**2
    
    # Verify by checking if d' = (dd/du) * u'
    # lhs is d' after substituting d = -u^2
    lhs = d_prime.subs(d, separatrix_candidate)
    
    # rhs is (d/du of the candidate) * u'
    d_du = sympy.diff(separatrix_candidate, u)
    rhs = d_du * u_prime
    
    # Check if LHS and RHS are equal
    if sympy.simplify(lhs - rhs) == 0:
        print("\nVerification successful:")
        print(f"The curve d = {separatrix_candidate} is an invariant manifold (a trajectory).")
        
        # Check if the curve passes through the saddle point
        d_saddle, u_saddle = saddle_point[d], saddle_point[u]
        if separatrix_candidate.subs(u, u_saddle) == d_saddle:
            print(f"The curve also passes through the saddle point ({d_saddle}, {u_saddle}).")
            print("Therefore, this curve is a separatrix for the system.\n")

            # Step 6: Print the final equation with explicit numbers
            print("--- Final Separatrix Equation ---")
            coeff = -1
            power = 2
            print(f"The equation of the separatrix is d(u) = {coeff} * u^{power}")
        else:
            print("The curve does not pass through the saddle point.")
    else:
        print("Verification failed: The curve is not a trajectory of the system.")


if __name__ == '__main__':
    find_separatrix()