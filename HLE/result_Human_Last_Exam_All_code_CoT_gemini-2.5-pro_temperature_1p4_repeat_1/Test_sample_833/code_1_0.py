# This script calculates the value of the expression for the specific stationary shock-front solution.
# The expression is E = (1-2u)*ubar*ubar_x, since u_t = 0 for a stationary solution.

# --- Case 1: Approaching the shock from the right (x -> 0+) ---
u_plus = 0.0
ubar_plus = -0.5
ubar_x_plus = 0.5

# Calculate E for x -> 0+
E_plus = (1 - 2 * u_plus) * ubar_plus * ubar_x_plus

print("Analysis for x -> 0+:")
print(f"u(0+) = {u_plus}")
print(f"ubar(0+) = {ubar_plus}")
print(f"ubar_x(0+) = {ubar_x_plus}")
print(f"The expression E(0+) evaluates to (1 - 2 * {u_plus}) * ({ubar_plus}) * ({ubar_x_plus}) = {E_plus}")
print("-" * 30)

# --- Case 2: Approaching the shock from the left (x -> 0-) ---
u_minus = 1.0
ubar_minus = -0.5
ubar_x_minus = -0.5

# Calculate E for x -> 0-
E_minus = (1 - 2 * u_minus) * ubar_minus * ubar_x_minus

print("Analysis for x -> 0-:")
print(f"u(0-) = {u_minus}")
print(f"ubar(0-) = {ubar_minus}")
print(f"ubar_x(0-) = {ubar_x_minus}")
print(f"The expression E(0-) evaluates to (1 - 2 * {u_minus}) * ({ubar_minus}) * ({ubar_x_minus}) = {E_minus}")
print("-" * 30)

# The infimum found using this special solution is the value calculated.
a = E_plus

print(f"The analysis of the shock-front solution shows that the expression can reach the value {a}.")
print(f"This indicates that the greatest lower bound is a = {a}.")