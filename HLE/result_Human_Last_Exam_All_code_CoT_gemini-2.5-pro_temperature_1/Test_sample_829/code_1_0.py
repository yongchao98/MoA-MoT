import numpy as np
from scipy.optimize import differential_evolution

# Define the component functions based on the variable u.
# h(u) corresponds to the local part of the flux F.
def h(u_val):
    return u_val * (1 - u_val)**2

# Derivatives of h(u) with respect to u.
def h_prime(u_val):
    # (u(1-u)^2)' = (u-2u^2+u^3)' = 1-4u+3u^2
    return 1 - 4 * u_val + 3 * u_val**2

def h_double_prime(u_val):
    # (1-4u+3u^2)' = -4+6u
    return -4 + 6 * u_val

def h_triple_prime(u_val):
    # (-4+6u)' = 6
    return 6.0

# Define the objective function to be maximized.
# The optimizer will minimize its negative.
# The input 'x' is a vector [u, u_bar, u_x1, u_bar_x1].
def calculate_expression(x):
    """Calculates the value of the expression E for a given state vector."""
    u, u_bar, u_x1, u_bar_x1 = x

    # Calculate the flux F and its value at x+1
    F_val = h(u) * np.exp(-u_bar)
    F_x1_val = h(u_x1) * np.exp(-u_bar_x1)

    # Calculate the derivatives of F with respect to u
    F_u = h_prime(u) * np.exp(-u_bar)
    F_uu = h_double_prime(u) * np.exp(-u_bar)
    # This is shorthand for d(F_uu)/du
    F_uuu = h_triple_prime(u) * np.exp(-u_bar)

    # The expression to be maximized, derived from the PDE
    # E = (u_x1 - u) * (F_uuu * F - F_uu * F_u) - F_uu * (F - F_x1)
    term1_factor = F_uuu * F_val - F_uu * F_u
    term1 = (u_x1 - u) * term1_factor
    term2 = F_uu * (F_val - F_x1_val)
    
    E = term1 - term2
    return E

def objective_function(x):
    """The function to be minimized (negative of the expression E)."""
    return -calculate_expression(x)

# Define the bounds for the four variables [u, u_bar, u_x1, u_bar_x1].
# Each variable is constrained to be in [0, 1].
bounds = [(0, 1), (0, 1), (0, 1), (0, 1)]

# Use differential_evolution for global optimization.
result = differential_evolution(objective_function, bounds, seed=42)

# The maximum value is the negative of the minimum value found by the optimizer.
max_value = -result.fun
# The state at which the maximum occurs.
optimal_point = result.x

# --- Output the results as requested ---
print("Finding the maximum of the expression E = (u_x+1 - u) * (F_uuu*F - F_uu*F_u) - F_uu * (F - F_x+1)\n")

u_opt, u_bar_opt, u_x1_opt, u_bar_x1_opt = optimal_point

# Recalculate components at the optimal point for printing
F_opt = h(u_opt) * np.exp(-u_bar_opt)
F_x1_opt = h(u_x1_opt) * np.exp(-u_bar_x1_opt)
F_u_opt = h_prime(u_opt) * np.exp(-u_bar_opt)
F_uu_opt = h_double_prime(u_opt) * np.exp(-u_bar_opt)
F_uuu_opt = h_triple_prime(u_opt) * np.exp(-u_bar_opt)

term1_factor_opt = F_uuu_opt * F_opt - F_uu_opt * F_u_opt

print(f"The maximum value is obtained near the point:")
print(f"u       = {u_opt:.4f}")
print(f"u_bar   = {u_bar_opt:.4f}")
print(f"u_x+1   = {u_x1_opt:.4f}")
print(f"u_bar_x+1 = {u_bar_x1_opt:.4f}\n")
# Analytically, the maximum is at u=0, u_bar=0, u_x+1=1. We round for clarity.
print("This corresponds to the analytical solution at u=0, u_bar=0, u_x+1=1.\n")


print("Plugging the analytical values into the equation:")
# Use analytical point for clean calculation breakdown
u_a, u_bar_a, u_x1_a = 0.0, 0.0, 1.0
F_a = h(u_a) * np.exp(-u_bar_a)
F_x1_a = h(u_x1_a) * np.exp(-0) # u_bar_x1 does not matter since h(1)=0
F_u_a = h_prime(u_a) * np.exp(-u_bar_a)
F_uu_a = h_double_prime(u_a) * np.exp(-u_bar_a)
F_uuu_a = h_triple_prime(u_a) * np.exp(-u_bar_a)
term1_factor_a = F_uuu_a * F_a - F_uu_a * F_u_a

print(f"F(u=0, u_bar=0) = {F_a:.4f}")
print(f"F(u_x+1=1) = {F_x1_a:.4f}")
print(f"F_u(u=0, u_bar=0) = {F_u_a:.4f}")
print(f"F_uu(u=0, u_bar=0) = {F_uu_a:.4f}")
print(f"F_uuu(u=0, u_bar=0) = {F_uuu_a:.4f}\n")

print("E = ({u_x1:.1f} - {u:.1f}) * ({F_uuu:.1f}*{F:.1f} - {F_uu:.1f}*{F_u:.1f}) - {F_uu:.1f} * ({F:.1f} - {F_x1:.1f})".format(
    u_x1=u_x1_a, u=u_a, F_uuu=F_uuu_a, F=F_a, F_uu=F_uu_a, F_u=F_u_a, F_x1=F_x1_a
))
print(f"E = ({u_x1_a - u_a:.1f}) * ({term1_factor_a:.1f}) - {F_uu_a:.1f} * ({F_a - F_x1_a:.1f})")
print(f"E = 1.0 * (6.0 * 0.0 - (-4.0) * 1.0) - (-4.0) * (0.0 - 0.0)")
print(f"E = 1.0 * 4.0 - 0.0")
print(f"E = {max_value:.4f}\n")

print("Final Answer:")
print(f"The maximum of the expression is {max_value:.1f}")
