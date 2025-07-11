import numpy as np
from scipy.optimize import fsolve

def solve_optimal_q(p_val):
    """
    Solves the first-order condition for the optimal search intensity q for a given p.
    The FOC is F(q, p) = 0, where:
    F(q, p) = exp(-2q) * [(1-p)(3-2q) + 2p(1-q)^2] - 2(1-q)
    """
    
    # Define the equation to be solved for q.
    # fsolve finds the root of this function.
    def foc_equation(q, p):
        # q must be a probability between 0 and 1.
        if q <= 0 or q >= 1:
            return 1e9  # Return a large value to guide the solver away from invalid regions.
        
        term_exp = np.exp(-2 * q)
        term_poly = (1 - p) * (3 - 2 * q) + 2 * p * (1 - q)**2
        lhs = term_exp * term_poly
        rhs = 2 * (1 - q)
        
        return lhs - rhs

    # An initial guess for q. Since it's a probability, 0.5 is a reasonable start.
    q_initial_guess = 0.5
    
    # Use fsolve to find the root of the equation.
    # The 'args' parameter passes the value of p to the equation.
    q_solution, = fsolve(foc_equation, q_initial_guess, args=(p_val,))
    
    return q_solution

# --- Main Calculation ---
# Problem parameters
p0 = 0.2
# A small value for h to calculate the numerical derivative
h = 1e-6

# 1. Calculate the optimal q at p = 0.2
q0 = solve_optimal_q(p0)

# 2. Calculate the numerical derivative dq/dp
q_plus_h = solve_optimal_q(p0 + h)
q_minus_h = solve_optimal_q(p0 - h)
dq_dp = (q_plus_h - q_minus_h) / (2 * h)

# 3. Calculate the elasticity E = (dq/dp) * (p/q)
elasticity = dq_dp * (p0 / q0)

# --- Output the results ---
# The problem asks to output the numbers in the final equation.
# The equation for elasticity is: Elasticity = (dq/dp) * (p / q)
print("Calculating the elasticity of optimal q with respect to p.")
print("The elasticity is calculated using the formula: E = (dq/dp) * (p / q)")
print("\nComponent values at p = 0.2:")
print(f"Optimal search intensity q = {q0:.5f}")
print(f"Probability of unemployment p = {p0:.5f}")
print(f"Derivative (dq/dp) = {dq_dp:.5f}")

# Final equation with numbers plugged in
print(f"\nElasticity = {dq_dp:.5f} * ({p0:.5f} / {q0:.5f})")
print(f"Final calculated elasticity = {elasticity:.5f}")
print(f"\nThe elasticity of optimal q with respect to p around p=0.2, rounded to three decimals, is: {elasticity:.3f}")

print(f"\n<<<{elasticity:.3f}>>>")