import math

# --- Plan ---
# 1. Define the physical constants for the problem.
# 2. Calculate the initial gauge pressure at the bottom of the river when the water is at rest.
# 3. Use Bernoulli's principle to set up the equation: The initial pressure energy is converted
#    to kinetic energy as the final pressure drops to zero.
#    P_initial = P_final + (1/2) * rho * v^2
# 4. Substitute the values and solve for the velocity v.
# 5. Print the equation with the numbers substituted in, and then print the final answer.

# --- Code Execution ---

# 1. Define constants
rho = 1000  # Density of water in kg/m^3
g = 10      # Acceleration due to gravity in m/s^2
H = 10      # Depth of the river in m

# 2. Calculate initial pressure at the bottom (when v=0)
P_initial = rho * g * H

# The final pressure at the bottom is given as 0
P_final = 0

# 3. Set up the Bernoulli equation. The initial pressure is converted to kinetic energy.
# P_initial = P_final + 0.5 * rho * v^2
# P_initial - P_final = 0.5 * rho * v^2
# We are solving for v.

# 4. Solve for v
v_squared = (P_initial - P_final) / (0.5 * rho)
v = math.sqrt(v_squared)

# 5. Print the equation and the final answer
print("From Bernoulli's principle, we relate the initial and final states:")
print("P_final + 1/2 * rho * v^2 = P_initial")
print("\nSubstituting the given values into the equation:")
# We explicitly show each number in the final equation.
print(f"{P_final} + 1/2 * {rho} * v^2 = {rho} * {g} * {H}")
print("\nSimplifying the equation:")
print(f"{P_final} + {0.5 * rho} * v^2 = {P_initial}")
print("\nSolving for v:")
print(f"v = sqrt(({P_initial} - {P_final}) / {0.5 * rho})")
print(f"The required flow speed is {v} m/s.")
