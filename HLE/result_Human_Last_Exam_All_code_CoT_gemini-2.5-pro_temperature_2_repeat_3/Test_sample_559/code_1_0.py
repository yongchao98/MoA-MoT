import sympy as sp

def solve_separatrix():
    """
    This function verifies the proposed equation for the separatrix of the given system of ODEs.
    
    The system is:
    d'(t) = 2d^2(t) + (-3u(t) + 5u^2(t))d(t) - u(t)(1-u(t))u^2(t)
    u'(t) = (u(t)-1)u^2(t)
    
    Through analysis, we find the equilibrium points at (0,0), (1,0), and (1,-1).
    The point (1,-1) is a saddle point. The separatrix is often the unstable
    manifold of a saddle point. We propose that the separatrix connecting the
    saddle point (1,-1) to the critical point (0,0) is given by the equation d = -u^2.
    
    This script will verify this claim symbolically.
    """
    
    # Define symbols and functions for symbolic manipulation
    t = sp.Symbol('t')
    u = sp.Function('u')(t)
    d = sp.Function('d')(t)

    # Define the original system of ODEs
    d_prime_rhs = 2*d**2 + (-3*u + 5*u**2)*d - u**3*(1-u)
    u_prime = u**3 - u**2

    # Define the proposed separatrix equation: d(t) = -u(t)^2
    proposed_d_solution = -u**2
    
    print("Verifying the proposed separatrix equation: d = -u^2\n")

    # 1. Calculate the LHS of the first ODE: d'(t)
    # We differentiate the proposed solution d(t) = -u(t)^2 with respect to t
    # using the chain rule: d/dt(-u^2) = -2u * du/dt
    lhs = sp.diff(proposed_d_solution, t)

    # Now substitute the expression for u'(t) from the second ODE
    lhs_final = lhs.subs(sp.Derivative(u, t), u_prime)
    
    # 2. Calculate the RHS of the first ODE
    # Substitute d = -u^2 into the expression for d'
    rhs_final = d_prime_rhs.subs(d, proposed_d_solution)

    # 3. Simplify both sides and check for equality
    lhs_simplified = sp.simplify(lhs_final)
    rhs_simplified = sp.simplify(rhs_final)

    print(f"Left-Hand Side (d'(t)) simplifies to: {lhs_simplified}")
    print(f"Right-Hand Side of the equation for d' simplifies to: {rhs_simplified}")

    # The verification consists of checking if LHS and RHS are identical.
    if sp.simplify(lhs_final - rhs_final) == 0:
        print("\nVerification successful: The Left-Hand Side is equal to the Right-Hand Side.")
        print("This confirms that d = -u^2 is an exact solution (trajectory) of the system.")
    else:
        print("\nVerification failed: The proposed solution is not correct.")

    # Print the final answer with its constituent numbers
    a = -1
    b = 2
    print("\nThe equation of the separatrix is d = a * u**b")
    print("The numbers in the final equation are:")
    print(f"Coefficient 'a': {a}")
    print(f"Exponent 'b': {b}")
    print(f"\nFinal Equation: d = {a}*u**{b}")

solve_separatrix()
<<<d = -u^2>>>