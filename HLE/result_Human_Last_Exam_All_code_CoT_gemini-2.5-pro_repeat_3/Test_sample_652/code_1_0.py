import numpy as np
from scipy.optimize import fsolve

def solve_for_q(p):
    """
    Solves the system of equations for a given p to find the optimal q.
    """
    # Equation 1: q = (ln(20) - ln(b)) / 2
    # This can be inverted to b = 20 * exp(-2q)
    
    # We are looking for the root of the following function of q:
    # f(q) = 20*exp(-2q) - 20*(1-q) / [ (1.5-q) + p*(1-q)^2/(1-p) ]
    def root_function(q, p_val):
        if not (0 < q < 1):
            return 1e9 # Penalize values outside the valid range
        
        b_from_q = 20 * np.exp(-2 * q)
        
        numerator = 20 * (1 - q)
        denominator = (1.5 - q) + p_val * (1 - q)**2 / (1 - p_val)
        
        if abs(denominator) < 1e-9:
             return 1e9

        b_from_govt = numerator / denominator
        
        return b_from_q - b_from_govt

    # Initial guess for q. If b=10, q = ln(2)/2 approx 0.35
    initial_guess_q = 0.35
    q_solution, = fsolve(root_function, initial_guess_q, args=(p,))
    return q_solution

# Define the parameters for the elasticity calculation
p = 0.2
dp = 1e-6 # A small change in p for the numerical derivative

# Calculate the optimal q for p and p+dp
q_star = solve_for_q(p)
q_new = solve_for_q(p + dp)

# Calculate the numerical derivative dq/dp
dq_dp = (q_new - q_star) / dp

# Calculate the elasticity
elasticity = dq_dp * (p / q_star)

# Output the final equation with the calculated numbers
print("The elasticity is calculated using the formula: ε = (dq/dp) * (p/q)")
print("The derivative dq/dp is approximated numerically as (q_new - q_star) / dp")
print("\n--- Calculation Steps ---")
print(f"For p = {p:.6f}, the optimal search intensity is q* = {q_star:.6f}")
print(f"For p_new = {p + dp:.6f}, the optimal search intensity is q_new = {q_new:.6f}")
print(f"The change in p is dp = {dp}")
print(f"The change in q is dq = {q_new - q_star:.6e}")
print(f"The derivative dq/dp ≈ ({q_new:.6f} - {q_star:.6f}) / {dp:.6f} = {dq_dp:.4f}")
print(f"The elasticity ε ≈ {dq_dp:.4f} * ({p} / {q_star:.6f}) = {elasticity:.3f}")
print("---")

<<<1.937>>>