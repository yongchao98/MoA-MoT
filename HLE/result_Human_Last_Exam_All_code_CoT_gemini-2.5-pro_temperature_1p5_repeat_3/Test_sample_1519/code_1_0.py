import math

# --- Given Parameters ---
# Density of water (rho) in kg/m^3
rho = 1000
# Acceleration due to gravity (g) in m/s^2
g = 10
# Depth of the river (H) in meters
H = 10

# --- Physics and Equation ---
# According to Bernoulli's principle, the pressure in the flowing water is reduced
# by the dynamic pressure (1/2 * rho * v^2).
# The pressure at the bottom when flowing is: P_flowing = P_static - P_dynamic
# P_static is the pressure at rest: rho * g * H
# We want to find the speed 'v' where P_flowing = 0.
# So, the equation is: 0 = (rho * g * H) - (1/2 * rho * v^2)

# Rearranging the equation to solve for v:
# 1/2 * rho * v^2 = rho * g * H

# --- Calculation ---
# We can solve for v from the rearranged equation.
# v^2 = 2 * g * H
v_squared = 2 * g * H
v = math.sqrt(v_squared)

# --- Output ---
print("The relationship between static pressure, dynamic pressure, and the resulting pressure is:")
print("Pressure_flowing = (rho * g * H) - (1/2 * rho * v^2)")
print("\nTo find the speed 'v' where the pressure at the bottom becomes 0, we solve the equation:")
print(f"0 = ({rho} * {g} * {H}) - (1/2 * {rho} * v^2)")
print("\nThis can be rearranged as:")
print(f"1/2 * {rho} * v^2 = {rho} * {g} * {H}")

# Calculate the values for the final equation line
p_static_value = rho * g * H
half_rho_value = 0.5 * rho
print(f"Substituting the values, we get:")
print(f"{half_rho_value} * v^2 = {p_static_value}")

print("\nSolving for v:")
print(f"v^2 = 2 * {g} * {H}")
print(f"v^2 = {v_squared}")
print(f"v = sqrt({v_squared})")
print(f"\nThe required flow speed is {v:.2f} m/s.")