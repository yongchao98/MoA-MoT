import numpy as np
from scipy.optimize import minimize

def get_q(b):
    """Calculates the worker's optimal re-employment probability q for a given benefit b."""
    # Ensure b is valid for log and results in q > 0
    if b <= 0 or b >= 20:
        return -np.inf # Invalid b
    return (np.log(20) - np.log(b)) / 2

def objective_function(b, p):
    """
    The government's objective function (worker's expected utility) to be maximized.
    We return its negative because we will use a minimizer.
    """
    # Calculate q based on the worker's optimization
    q = get_q(b)
    
    # Check if q is within the valid range (0, 1)
    # The optimizer bounds should handle this, but it's good practice for robustness
    if not (0 < q < 1):
        return np.inf

    # Calculate tax t from the government's budget constraint: (1-p)t = p(1-q)b
    if p == 0 or p >= 1:
        # Handle edge cases for p
        return np.inf
    t = (p * b * (1 - q)) / (1 - p)
    
    # Consumption in the initially employed state must be positive
    consumption_employed = 20 - t
    if consumption_employed <= 0:
        return np.inf

    # Worker's maximized expected utility in the unemployed state is ln(b) + q^2
    # This simplification comes from substituting the worker's FOC for q back into their utility.
    utility_unemployed_state = np.log(b) + q**2
    
    # Government's objective: total expected utility for the worker
    total_expected_utility = (1 - p) * np.log(consumption_employed) + p * utility_unemployed_state
    
    # We use a minimizer, so we return the negative of the utility
    return -total_expected_utility

def find_optimal_q(p):
    """
    Finds the optimal q by first finding the optimal b that maximizes
    the worker's expected utility for a given p.
    """
    # Initial guess for b. A value within the bounds is a good starting point.
    initial_b = 10.0
    
    # Bounds for b to ensure q is in (0, 1)
    # q > 0 => b < 20
    # q < 1 => b > 20/e^2
    b_bounds = [(20 / np.exp(2) + 1e-9, 20 - 1e-9)]
    
    # Minimize the negative expected utility to find the optimal b
    result = minimize(
        fun=objective_function,
        x0=initial_b,
        args=(p,),
        bounds=b_bounds,
        method='L-BFGS-B'
    )
    
    if not result.success:
        raise ValueError(f"Optimization failed for p={p}: {result.message}")
        
    optimal_b = result.x[0]
    
    # Calculate the corresponding optimal q
    optimal_q = get_q(optimal_b)
    
    return optimal_q

if __name__ == '__main__':
    # Given probability of unemployment
    p = 0.2
    
    # Small change in p for numerical differentiation
    dp = 1e-6

    # 1. Calculate the optimal q for the given p
    q_at_p = find_optimal_q(p)
    
    # 2. Calculate the optimal q for p + dp
    q_at_p_plus_dp = find_optimal_q(p + dp)
    
    # 3. Calculate the numerical derivative dq/dp
    dq_dp = (q_at_p_plus_dp - q_at_p) / dp
    
    # 4. Calculate the elasticity of q with respect to p
    elasticity = dq_dp * (p / q_at_p)

    # Print the numbers in the final equation
    print("Elasticity of optimal q with respect to p:")
    print(f"The calculation is: (dq/dp) * (p / q)")
    print(f"p = {p}")
    print(f"Optimal q at p={p}: {q_at_p:.6f}")
    print(f"Derivative dq/dp at p={p}: {dq_dp:.6f}")
    print(f"Elasticity = {dq_dp:.6f} * ({p} / {q_at_p:.6f})")
    print(f"Final Answer (3 decimals): {elasticity:.3f}")
    print(f"\n<<< {elasticity:.3f} >>>")