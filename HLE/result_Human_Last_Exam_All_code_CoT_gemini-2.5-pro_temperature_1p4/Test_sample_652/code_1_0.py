import numpy as np
from scipy.optimize import brentq

def solve_and_calculate_elasticity():
    """
    This function solves the economic model to find the elasticity of job search
    intensity (q) with respect to the probability of unemployment (p).
    """
    # Define model parameters
    w = 20.0  # Wage
    p = 0.2   # Probability of becoming unemployed

    # --- Step 1: Numerically solve for b and q at p=0.2 ---

    # We need to solve the non-linear equation derived from the government's FOC:
    # b * (1.5 - q - 0.5*p) = 20
    # where q is also a function of b: q = (ln(20) - ln(b)) / 2
    # We define a function f(b) whose root we want to find.
    def f(b, p_val, w_val):
        # q is constrained to be in (0, 1), which implies b must be in (w/e^2, w).
        # We handle cases outside this range to guide the solver.
        if b <= 0:
            return -np.inf
        # ln(b) is only defined for b > 0
        q = (np.log(w_val) - np.log(b)) / 2
        # The FOC equation: b * (1.5 - q - 0.5*p) - w = 0
        return b * (1.5 - q - 0.5 * p_val) - w_val

    # The valid range for b is (w/e^2, w), which is approx (2.7, 20).
    # We search for the root within this interval.
    lower_bound = w / np.exp(2) + 1e-9
    upper_bound = w - 1e-9
    
    # Use the Brent method to find the root b_val
    try:
        b_val = brentq(f, lower_bound, upper_bound, args=(p, w))
    except ValueError:
        print("Could not find a solution for b in the given interval.")
        return

    # Calculate the corresponding optimal q
    q_val = (np.log(w) - np.log(b_val)) / 2

    # --- Step 2: Calculate the elasticity ---

    # The formula for dq/dp was derived using implicit differentiation:
    # dq/dp = -0.25 * b / (w + 0.5 * b)
    dq_dp = -0.25 * b_val / (w + 0.5 * b_val)

    # The elasticity E = (dq/dp) * (p/q)
    elasticity = dq_dp * (p / q_val)

    # --- Step 3: Print the results as requested ---
    
    print("This script calculates the elasticity of optimal job search intensity (q) with respect to the probability of unemployment (p).")
    print(f"\nGiven parameters: p = {p}, w = {w}")
    print("\nIntermediate values at the optimum:")
    print(f"Benefit (b) = {b_val:.4f}")
    print(f"Search Intensity (q) = {q_val:.4f}")
    
    print("\nFinal elasticity calculation:")
    # "output each number in the final equation!"
    print(f"E = (dq/dp) * (p/q)")
    print(f"E = ({dq_dp:.4f}) * ({p:.1f} / {q_val:.4f})")
    print(f"E = {elasticity:.3f}")

# Execute the function
solve_and_calculate_elasticity()