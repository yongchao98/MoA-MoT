import sympy as sp

def find_separatrix():
    """
    This function demonstrates the process of finding and verifying the separatrix
    for the given system of differential equations.
    """
    # 1. Define symbolic variables for the verification.
    # We are looking for a relationship between d and u, i.e., d=f(u).
    u, d = sp.symbols('u d')

    # 2. Define the expressions for the derivatives from the ODE system.
    # d'(t) = 2*d**2 + (-3*u + 5*u**2)*d - u**3*(1-u)
    # u'(t) = u**2*(u-1)
    d_prime_ode_expr = 2*d**2 + (-3*u + 5*u**2)*d - u**3*(1-u)
    u_prime_ode_expr = u**2 * (u - 1)

    # 3. State the proposed separatrix equation, which is found through analysis.
    # The analysis (finding equilibria, identifying the saddle at (-1, 1),
    # and testing for an invariant manifold) leads to the hypothesis d = -u^2.
    d_solution = -u**2

    print("Step 1: State the problem and the proposed solution.")
    print("System of ODEs:")
    print(f"d'(t) = 2*d**2 + (-3*u + 5*u**2)*d - u*(-u + 1)*u**2")
    print(f"u'(t) = (u - 1)*u**2")
    print("\nThrough analysis, we find that (-1, 1) is a saddle point.")
    print(f"We hypothesize the separatrix connected to it is the curve d = {d_solution}.")
    print("\nStep 2: Verify the proposed solution.")
    print("The verification checks if a trajectory on the curve d = -u^2 satisfies the system.")
    print("This is true if d'(t) calculated from the chain rule matches d'(t) from the first ODE.")

    # 4. Verification calculation.
    # Method A: Calculate d'(t) from the chain rule: d'(t) = (d(d)/du) * u'(t).
    # The derivative of d_solution with respect to u is sp.diff(d_solution, u).
    d_prime_from_chain_rule = sp.diff(d_solution, u) * u_prime_ode_expr
    d_prime_from_chain_rule = sp.expand(d_prime_from_chain_rule)
    
    print("\n   Calculation from Chain Rule:")
    print(f"   d'(t) = d/dt(-u**2) = (d/du(-u**2)) * u'(t) = ({sp.diff(d_solution, u)}) * ({u_prime_ode_expr})")
    print(f"   d'(t) simplified = {d_prime_from_chain_rule}")

    # Method B: Calculate d'(t) from the d-ODE by substituting d = d_solution.
    d_prime_from_ode = d_prime_ode_expr.subs(d, d_solution)
    d_prime_from_ode = sp.expand(d_prime_from_ode)

    print("\n   Calculation from the d-ODE:")
    print(f"   d'(t) = 2*({d_solution})**2 + (-3*u + 5*u**2)*({d_solution}) - u**3*(1-u)")
    print(f"   d'(t) simplified = {d_prime_from_ode}")
    
    # 5. Conclude the verification.
    print("\nStep 3: Conclusion of Verification.")
    if sp.simplify(d_prime_from_chain_rule - d_prime_from_ode) == 0:
        print("Result: The two expressions for d'(t) are identical.")
        print("This confirms that d = -u^2 is a solution trajectory and therefore a separatrix of the system.\n")
    else:
        print("Result: Verification failed. The expressions are not identical.\n")

    # 6. Print the final equation as requested.
    print("Final Answer: The Equation of the Separatrix")
    print("The equation is d = -u^2, or equivalently, d + u^2 = 0.")
    
    # Define and print each number in the final equation d + c*u^p = C
    c = 1
    p = 2
    C = 0
    print(f"This can be written in the form d + c * u^p = C, where:")
    print(f"d + ({c}) * u**({p}) = {C}")

if __name__ == '__main__':
    find_separatrix()