import sympy

def find_separatrix():
    """
    Finds and verifies the separatrix for the given system of differential equations.
    """
    # Define symbolic variables
    d, u = sympy.symbols('d u')

    # Define the system of ODEs
    # d'(t) = 2d^2 + (-3u + 5u^2)d - u^3(1-u)
    # u'(t) = u^2(u-1)
    d_prime = 2*d**2 + (-3*u + 5*u**2)*d - u**3*(1-u)
    u_prime = u**2 * (u - 1)

    print("System of ODEs:")
    print(f"d'/dt = {d_prime}")
    print(f"u'/dt = {u_prime}\n")

    # --- Step 1 & 2: Find and classify equilibrium points ---
    print("--- Finding and Classifying Equilibrium Points (for u>=0) ---")
    
    # Solve for equilibrium points where d'=0 and u'=0
    equilibria = sympy.solve([d_prime, u_prime], (d, u), dict=True)
    
    saddle_point = None
    for eq in equilibria:
        # Per problem statement, we are on the u >= 0 half-plane
        if eq[u] >= 0:
            print(f"Found equilibrium point: (d, u) = ({eq[d]}, {eq[u]})")
            
            # Linearize the system by computing the Jacobian matrix
            F = sympy.Matrix([d_prime, u_prime])
            var_vec = sympy.Matrix([d, u])
            J = F.jacobian(var_vec)
            
            # Evaluate the Jacobian at the equilibrium point
            J_eq = J.subs(eq)
            eigenvals = J_eq.eigenvals()

            # Classify the point based on eigenvalues
            print(f"  Eigenvalues at this point: {list(eigenvals.keys())}")
            if J_eq.is_zero:
                print("  Classification: Non-hyperbolic (degenerate)")
            else:
                real_parts = [sympy.re(ev) for ev in eigenvals.keys()]
                if any(r > 0 for r in real_parts) and any(r < 0 for r in real_parts):
                    saddle_point = eq
                    print("  Classification: Saddle Point")
                elif all(r > 0 for r in real_parts):
                    print("  Classification: Unstable Node (Source)")
                else: # Covers stable nodes, spirals, etc.
                     print("  Classification: Other/Stable")


    if not saddle_point:
        print("\nNo saddle point found. This method to find the separatrix is not applicable.")
        return

    print(f"\nA saddle point is identified at (d, u) = ({saddle_point[d]}, {saddle_point[u]}).\n")
    
    # --- Step 3: Propose and verify the separatrix equation ---
    print("--- Finding the Separatrix Equation ---")
    print("A separatrix is an invariant curve passing through a saddle point.")
    print("Let's test if the curve d = -u^2 is the separatrix.")

    # The proposed separatrix equation
    proposed_d = -u**2
    
    # To verify, we must show that if d(t) and u(t) are on this curve,
    # their derivatives d'(t) and u'(t) satisfy the system equations.
    # We substitute d = -u^2 into the first ODE and check for consistency.
    # LHS is d'(t). Using the chain rule: d'(t) = (d(d)/du) * u'(t).
    dd_du = sympy.diff(proposed_d, u)
    lhs_d_prime = dd_du * u_prime
    
    # RHS is the original expression for d', with d = -u^2 substituted.
    rhs_d_prime = d_prime.subs(d, proposed_d)
    
    print(f"\nVerifying if d = -u^2 is an invariant curve...")
    print(f"  d'(t) derived from the curve d={proposed_d} is: {sympy.simplify(lhs_d_prime)}")
    print(f"  d'(t) from the system's equation is: {sympy.simplify(rhs_d_prime)}")

    if sympy.simplify(lhs_d_prime - rhs_d_prime) == 0:
        print("\nVerification successful! The curve d = -u^2 is an invariant manifold of the system.")
        
        # Check if it passes through the saddle point
        saddle_d_val = proposed_d.subs(u, saddle_point[u])
        if saddle_d_val == saddle_point[d]:
             print(f"The curve also passes through the saddle point ({saddle_point[d]}, {saddle_point[u]}).")
             print("Therefore, it is a separatrix.")
             
             # Final Answer Formatting as requested
             print(f"\nThe equation for the separatrix is d = -u^2.")
             print(f"To explicitly show each number in the final equation: d = (-1) * u**2")
        else:
            print("However, the curve does not pass through the identified saddle point.")
    else:
        print("\nVerification failed. d = -u^2 is not the separatrix.")

# Execute the analysis
find_separatrix()