import math

# Define the given physical constants and initial conditions
p_static = 100000  # Initial static pressure at the bottom in N/m^2
rho = 1000         # Density of water in kg/m^3

# According to the problem, the pressure in the flowing river is P_flowing = P_static - (1/2)*rho*v^2
# We need to find the speed v for which P_flowing = 0.
# The equation to solve is:
# 0 = P_static - (1/2) * rho * v^2
# which rearranges to:
# P_static = (1/2) * rho * v^2

# Now, we solve for v.
# v^2 = 2 * P_static / rho
v_squared = (2 * p_static) / rho

# v = sqrt(v^2)
v = math.sqrt(v_squared)

# Print the final equation with all the numerical values plugged in.
print("The final equation we solve is P_static = (1/2) * rho * v^2")
print(f"By plugging in the values, we get:")
print(f"{p_static} = (1/2) * {rho} * v^2")
print("\nSolving this for v yields:")
print(f"v = sqrt(2 * {p_static} / {rho})")
print(f"v = {v:.2f} m/s")

# The exact answer is sqrt(200)
# print(f"<<<{math.sqrt(200)}>>>")