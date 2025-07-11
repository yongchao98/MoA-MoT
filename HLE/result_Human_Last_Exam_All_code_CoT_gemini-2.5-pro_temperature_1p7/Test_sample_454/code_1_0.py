import math

# --- 1. Define Parameters ---
# The problem gives parameters relative to k_B*T, so we can work with these ratios.
mu_over_kT = 0.1
epsilon_over_kT = -1 / (2 * math.pi)

# Coordination numbers
z_horizontial = 4
z_vertical = 8
# The total coordination number z for a site in a multilayer system (not on the surface layer)
# is the sum of its horizontal and vertical nearest neighbors.
z = z_horizontial + z_vertical

# --- 2. Solve the Self-Consistency Equation Numerically ---
# The equation is: n = 1 / (exp((z * epsilon/kT) * n - mu/kT) + 1)
# We use fixed-point iteration to find the value of n (average occupancy).

# Initial guess for the average occupancy n
n = 0.5

# Iterate to find the stable solution. 100 iterations are sufficient for convergence.
for _ in range(100):
    exponent = (z * epsilon_over_kT * n) - mu_over_kT
    n = 1.0 / (math.exp(exponent) + 1)

# Round the final answer to three decimal places as requested
final_n = round(n, 3)

# --- 3. Display the Final Equation and Result ---
# The prompt requires displaying the final equation with all numbers plugged in.
print("The self-consistency equation for the average occupancy <n> is:")
print("<n> = 1 / (exp((z * ε/k_B T) * <n> - μ/k_B T) + 1)")
print("\nAfter substituting the parameter values:")
print(f"z = {z}")
print(f"ε/k_B T ≈ {epsilon_over_kT:.3f}")
print(f"μ/k_B T = {mu_over_kT}")

print("\nThe final equation with the solution plugged in is:")

# We use the unrounded value of epsilon_over_kT for better precision in the intermediate step
term1 = z * epsilon_over_kT
term2 = term1 * final_n
full_exponent = term2 - mu_over_kT

# Final equation printout as requested
print(f"{final_n} = 1 / (exp(({z} * {epsilon_over_kT:.3f}) * {final_n} - {mu_over_kT}) + 1)")
print(f"{final_n} = 1 / (exp({term1:.3f} * {final_n} - {mu_over_kT}) + 1)")
print(f"{final_n} = 1 / (exp({term2:.3f} - {mu_over_kT}) + 1)")
print(f"{final_n} = 1 / (exp({full_exponent:.3f}) + 1)")

print(f"\nThe average occupancy per site <n> is {final_n}")
print("<<<0.848>>>")