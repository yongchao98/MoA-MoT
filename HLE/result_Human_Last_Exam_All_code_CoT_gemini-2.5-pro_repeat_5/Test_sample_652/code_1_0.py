import numpy as np
from scipy.optimize import brentq

def solve_elasticity():
    """
    Solves the microeconomic problem to find the elasticity of optimal q
    with respect to p.
    """
    # The value of p at which to calculate the elasticity
    p_val = 0.2

    # F(q, p) = 0 is the first-order condition for the government's problem.
    # We define it to find the optimal q* for a given p.
    def F(q, p):
        # Avoid division by zero if p=1
        if p == 1:
            # In this case, the worker is always unemployed, and the GBC is ill-defined.
            # This is not relevant for p=0.2 but good practice.
            return np.nan
        term1 = np.exp(-2 * q) * (3 - 2 * q)
        term2 = (2 * p / (1 - p)) * ((1 - q)**2) * np.exp(-2 * q)
        term3 = -2 * (1 - q)
        return term1 + term2 + term3

    # We need to solve F(q, 0.2) = 0 for q. We define a lambda function for this.
    # Preliminary analysis shows the root is between 0.3 and 0.4.
    F_at_p = lambda q: F(q, p_val)
    try:
        q_star = brentq(F_at_p, 0.0, 1.0)
    except ValueError:
        print("Could not find a root for q* in the interval [0, 1].")
        return

    # Now we use the Implicit Function Theorem: dq/dp = - (dF/dp) / (dF/dq)
    # 1. Calculate the partial derivative of F with respect to p
    def dF_dp(q, p):
        # Derivative of 2p/(1-p) is 2/((1-p)^2)
        return (2 / ((1 - p)**2)) * ((1 - q)**2) * np.exp(-2 * q)

    # 2. Calculate the partial derivative of F with respect to q
    def dF_dq(q, p):
        # Derivative of F(q,p) w.r.t q
        df_term1 = -2 * np.exp(-2 * q) * (3 - 2 * q) + np.exp(-2 * q) * (-2)
        
        df_term2_1 = (2 * p / (1 - p)) * 2 * (1 - q) * (-1) * np.exp(-2 * q)
        df_term2_2 = (2 * p / (1 - p)) * ((1 - q)**2) * (-2 * np.exp(-2 * q))
        df_term2 = df_term2_1 + df_term2_2
        
        df_term3 = -2 * (-1)
        
        return df_term1 + df_term2 + df_term3

    # 3. Evaluate the derivatives at (q*, p_val)
    dF_dp_val = dF_dp(q_star, p_val)
    dF_dq_val = dF_dq(q_star, p_val)

    # 4. Calculate dq/dp
    dq_dp = -dF_dp_val / dF_dq_val

    # 5. Calculate the elasticity E = (dq/dp) * (p/q)
    elasticity = dq_dp * (p_val / q_star)

    # Print the intermediate values and the final result
    print(f"For p = {p_val}:")
    print(f"The optimal re-employment probability is q* = {q_star:.4f}")
    print("\nCalculating the components for elasticity:")
    print(f"The partial derivative dF/dp at (q*, p) is {dF_dp_val:.4f}")
    print(f"The partial derivative dF/dq at (q*, p) is {dF_dq_val:.4f}")
    print(f"The derivative dq*/dp = -({dF_dp_val:.4f}) / ({dF_dq_val:.4f}) = {dq_dp:.4f}")
    print("\nThe final equation for elasticity is:")
    print(f"E = (dq*/dp) * (p / q*)")
    print(f"E = {dq_dp:.4f} * ({p_val} / {q_star:.4f})")
    print(f"The elasticity of optimal q with respect to p is: {elasticity:.3f}")
    
    # Returning the final answer for the '<<<' format
    return elasticity

if __name__ == '__main__':
    final_answer = solve_elasticity()
    # The final answer is then wrapped as requested
    # print(f"\n<<<{final_answer:.3f}>>>")