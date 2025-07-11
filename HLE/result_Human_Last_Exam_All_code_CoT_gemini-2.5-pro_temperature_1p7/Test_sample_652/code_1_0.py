import numpy as np
from scipy.optimize import brentq

def get_optimal_q(p):
    """
    Finds the optimal q for a given p by solving the first-order condition
    of the government's optimization problem.
    The FOC is:
    (exp(-2q) * (3 - 2q)) / (1 - (p/(1-p))*(1-q)*exp(-2q)) + 2(q-1) = 0
    """

    # Define the function whose root we want to find
    def foc_equation(q):
        # Avoid division by zero if p=1 or q is close to 1 causing issues
        if q >= 1.0:
            # The left part is positive, the right is 0, so overall > 0.
            # Return a large positive number to guide the solver.
            return 1e9

        ratio = p / (1 - p)
        exp_val = np.exp(-2 * q)
        
        numerator = exp_val * (3 - 2 * q)
        denominator = 1 - ratio * (1 - q) * exp_val
        
        # Check for non-positive denominator, which would imply c<0, invalid
        if denominator <= 0:
            # This combination of (p, q) is infeasible.
            # Return a large number to push the solver away.
            return 1e9
            
        return numerator / denominator + 2 * (q - 1)

    # Use a numerical solver (Brent's method) to find the root.
    # The bracket [0.001, 0.999] is a safe range for q.
    try:
        optimal_q = brentq(foc_equation, 0.001, 0.999, xtol=1e-12, rtol=1e-12)
        return optimal_q
    except ValueError:
        print(f"Error: Could not find a root for p={p} in the given interval.")
        return None

def calculate_elasticity():
    """
    Calculates the elasticity of optimal q with respect to p around p=0.2.
    """
    p1 = 0.2
    # Use a small change in p for numerical differentiation
    delta_p = 1e-7
    p2 = p1 + delta_p

    # Calculate optimal q for p1 and p2
    q1 = get_optimal_q(p1)
    q2 = get_optimal_q(p2)

    if q1 is None or q2 is None:
        return

    # Calculate the numerical derivative dq/dp
    dq_dp = (q2 - q1) / delta_p

    # Calculate the elasticity E = (dq/dp) * (p/q)
    elasticity = dq_dp * (p1 / q1)

    print(f"Calculating the elasticity E = (dq/dp) * (p/q):")
    print(f"Base probability p = {p1}")
    print(f"Optimal search intensity q at p={p1:.3f} is q1 = {q1:.6f}")
    print(f"Optimal search intensity q at p={p2:.7f} is q2 = {q2:.6f}")
    print(f"Numerical derivative dq/dp ≈ (q2 - q1) / delta_p = {dq_dp:.6f}")
    print(f"Elasticity = {dq_dp:.6f} * ({p1:.3f} / {q1:.6f}) = {elasticity:.3f}")
    
    # Final answer block
    print("\n<<<" + f"{elasticity:.3f}" + ">>>")


# Run the calculation
calculate_elasticity()
