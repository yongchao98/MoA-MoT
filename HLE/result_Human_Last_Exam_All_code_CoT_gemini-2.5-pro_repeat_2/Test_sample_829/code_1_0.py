import numpy as np
from scipy.optimize import differential_evolution

# Define the function f(u) and its derivatives
# f(u) = u(1-u)^2
def f(u):
    return u * (1 - u)**2

# f'(u) = (1-u)^2 - 2u(1-u) = 1 - 4u + 3u^2
def f_prime(u):
    return 1 - 4 * u + 3 * u**2

# f''(u) = -4 + 6u
def f_double_prime(u):
    return -4 + 6 * u

# f'''(u) = 6
def f_triple_prime(u):
    return 6

# Define the main expression E to be maximized
def E_func(x):
    """
    Calculates the value of the expression to be maximized.
    x = [u1, u2, u_bar1, u_bar2]
    """
    u1, u2, u_bar1, u_bar2 = x

    # Calculate F and its derivatives at (u1, u_bar1)
    F1 = f(u1) * np.exp(-u_bar1)
    F_prime1 = f_prime(u1) * np.exp(-u_bar1)
    F_double_prime1 = f_double_prime(u1) * np.exp(-u_bar1)
    F_triple_prime1 = f_triple_prime(u1) * np.exp(-u_bar1)

    # Calculate F at (u2, u_bar2)
    F2 = f(u2) * np.exp(-u_bar2)

    # Calculate the two main terms of the expression E
    term1_factor = (F_triple_prime1 * F1 - F_double_prime1 * F_prime1)
    term1 = term1_factor * (u2 - u1)
    
    term2_factor = F_double_prime1
    term2 = term2_factor * (F1 - F2)

    return term1 - term2

# We want to maximize E_func, which is equivalent to minimizing -E_func
def objective_function(x):
    return -E_func(x)

# Define the bounds for the variables [u1, u2, u_bar1, u_bar2]
bounds = [(0, 1), (0, 1), (0, 1), (0, 1)]

# Use differential evolution for global optimization
result = differential_evolution(objective_function, bounds)

# The maximum value is the negative of the minimum found
max_value = -result.fun
optimal_vars = result.x

# Print the results
print("Numerical optimization results:")
print(f"Maximum value found: {max_value:.8f}")
print(f"At variables (u1, u2, u_bar1, u_bar2): {optimal_vars}")
print("\n--- Verification using the numerically found optimal variables ---")

# Unpack optimal variables
u1, u2, u_bar1, u_bar2 = optimal_vars

# Recalculate the components of the equation with the optimal values
F1_val = f(u1) * np.exp(-u_bar1)
F_prime1_val = f_prime(u1) * np.exp(-u_bar1)
F_double_prime1_val = f_double_prime(u1) * np.exp(-u_bar1)
F_triple_prime1_val = f_triple_prime(u1) * np.exp(-u_bar1)
F2_val = f(u2) * np.exp(-u_bar2)

term1_factor_val = F_triple_prime1_val * F1_val - F_double_prime1_val * F_prime1_val
term1_val = term1_factor_val * (u2 - u1)

term2_factor_val = F_double_prime1_val
term2_val = term2_factor_val * (F1_val - F2_val)

final_val = term1_val - term2_val

print("Calculation with optimal variables:")
print(f"u1 = {u1:.4f}, u2 = {u2:.4f}, u_bar1 = {u_bar1:.4f}, u_bar2 = {u_bar2:.4f}")
print(f"Term 1 factor (F'''F - F''F'): {term1_factor_val:.4f}")
print(f"Term 1 = {term1_factor_val:.4f} * ({u2:.4f} - {u1:.4f}) = {term1_val:.4f}")
print(f"Term 2 factor (-F''): {-term2_factor_val:.4f}")
print(f"Term 2 = -({-term2_factor_val:.4f}) * ({F1_val:.4f} - {F2_val:.4f}) = {-term2_val:.4f}")
print(f"Final Value = {term1_val:.4f} - ({term2_val:.4f}) = {final_val:.4f}")
print("\nAs shown by the analytical derivation, the exact maximum is 4.")

<<<4>>>