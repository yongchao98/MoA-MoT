import numpy as np
from scipy.optimize import brentq

def solve_for_q(p):
    """
    Solves for the optimal q for a given probability p.

    The function to find the root of is derived from the government's FOC:
    (1-q)(20-t) = b(1.5 - q)
    where:
    b = 20 * exp(-2q)
    t = (p / (1-p)) * b * (1-q)
    """
    
    # Define the implicit function f(q) = 0 that determines the optimal q.
    def implicit_function(q):
        # The equation is valid for 0 < q < 1 and 0 < p < 1
        if q <= 0 or q >= 1 or p <= 0 or p >= 1:
            return float('nan')
        
        # Calculate b and t based on q and p
        b = 20 * np.exp(-2 * q)
        t = (p / (1 - p)) * b * (1 - q)

        # The equation from the first-order condition
        # (1-q)(20-t) - b(1.5-q) = 0
        return (1 - q) * (20 - t) - b * (1.5 - q)

    # Use a numerical solver (brentq) to find the root q.
    # We found through analysis that the root lies between 0.1 and 0.6.
    try:
        q_optimal = brentq(implicit_function, 0.01, 0.99)
        return q_optimal
    except ValueError:
        return float('nan')

def main():
    """
    Calculates the elasticity of optimal q with respect to p around p=0.2.
    """
    # Initial probability
    p0 = 0.2
    
    # Small change in p for numerical differentiation
    h = 1e-6
    p1 = p0 + h

    # Calculate optimal q for p0 and p1
    q0 = solve_for_q(p0)
    q1 = solve_for_q(p1)

    # Calculate the derivative dq/dp using finite difference
    dq_dp = (q1 - q0) / h

    # Calculate the elasticity
    elasticity = (dq_dp * p0) / q0
    
    # Print the results as requested
    print(f"The optimal re-employment probability (q) for p = {p0} is {q0:.5f}")
    print(f"The derivative of q with respect to p (dq/dp) is approximately {dq_dp:.5f}")
    print("\nThe elasticity of q with respect to p is calculated as:")
    print(f"ε = (dq/dp) * (p/q) = {dq_dp:.5f} * ({p0:.1f} / {q0:.5f}) = {elasticity:.3f}")
    
    # Final answer in the specified format
    print(f"\n<<<{elasticity:.3f}>>>")

if __name__ == "__main__":
    main()
