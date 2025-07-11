import numpy as np
from scipy.optimize import brentq

# Define the parameters of the model based on the problem description.
p = 0.2
w = 20

# The first-order condition (FOC) from the government's optimization problem, which defines
# the optimal job search intensity q* as a function of the unemployment probability p, is:
# exp(-2*q) * (3 - p - 2q) - 2 = 0
# For our given p = 0.2, this simplifies to:
# exp(-2*q) * (2.8 - 2*q) - 2 = 0
# We define this as a function to find its root, which will be our q*.
def find_q_func(q, p_val=p):
    """This is the FOC for the optimal q. We need to find the root of this function."""
    return np.exp(-2 * q) * (3 - p_val - 2 * q) - 2

# We solve for q* numerically in the interval (0, 1), as q is a probability.
# A quick check shows a unique root exists in this interval:
# At q=0, the function is 3 - p - 2 = 1 - p = 0.8 > 0
# As q->inf, the function approaches -2 < 0.
try:
    # Use the Brent method for robust root-finding.
    q_star = brentq(find_q_func, 0, 1)
except (ValueError, RuntimeError):
    q_star = None

if q_star is not None:
    # Using implicit differentiation on the FOC, we find the derivative of q* with respect to p:
    # dq*/dp = -exp(-2q*) / (4 + 2*exp(-2q*))
    exp_term = np.exp(-2 * q_star)
    dq_dp = -exp_term / (4 + 2 * exp_term)

    # Now, calculate the elasticity using its definition:
    # Elasticity E = (dq*/dp) * (p / q*)
    elasticity = dq_dp * (p / q_star)

    # Print the results, showing the components of the final calculation as requested.
    print("This script calculates the elasticity of the optimal job search intensity (q) with respect to the probability of unemployment (p).")
    print("\n--- Intermediate Values ---")
    print(f"Given probability of unemployment, p = {p}")
    print(f"Calculated optimal job search intensity, q* = {q_star:.5f}")
    print(f"Calculated derivative of q* with respect to p, dq*/dp = {dq_dp:.5f}")

    print("\n--- Final Elasticity Calculation ---")
    print(f"Elasticity E = (dq*/dp) * (p / q*)")
    print(f"E = {dq_dp:.5f} * ({p} / {q_star:.5f})")
    print(f"E = {elasticity:.5f}")

    # Final answer formatted as requested, rounded to three decimals.
    print(f"\nThe elasticity of optimal q with respect to p, rounded to three decimals, is {elasticity:.3f}.")
    print(f'<<<{elasticity:.3f}>>>')

else:
    print("Could not find the optimal q*. Check the function and interval.")
