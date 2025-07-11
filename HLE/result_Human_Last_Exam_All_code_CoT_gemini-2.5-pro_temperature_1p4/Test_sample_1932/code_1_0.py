import math

# Description of the reasoning
# The change in weight of the hourglass, Delta_W, can be estimated by considering the largest possible dynamic effect.
# This occurs at the beginning of the sand's flow.
# At this moment, the downward force from the impact of the first falling grains is maximal,
# while the counteracting effect from the weight of sand in the air is still negligible.
# The impact force is given by the mass flow rate (m_dot) multiplied by the impact velocity (v_f).
# m_dot = total_mass / total_time = (rho * V_sand) / t = (rho * (pi * d**2 / 4) * h) / t
# The largest impact velocity occurs for the largest fall height, H.
# v_f = sqrt(2 * g * H)
# Thus, the largest weight change is Delta_W = m_dot * v_f.
# Delta_W = (pi * d**2 * h * rho) / (4 * t) * sqrt(2 * g * H)

# Given parameters
d_cm = 1.0  # diameter in cm
h_cm = 2.0  # sand column height in cm
H_cm = 4.0  # chamber height in cm
rho = 1500  # density in kg/m^3
t_min = 1.0  # time in minutes
g = 9.8  # acceleration due to gravity in m/s^2

# Convert units to SI
d = d_cm / 100.0  # m
h = h_cm / 100.0  # m
H = H_cm / 100.0  # m
t = t_min * 60.0  # s

# Chosen formula
formula_str = "Delta_W = (pi * d**2 * h * rho / (4 * t)) * sqrt(2 * g * H)"

# Print the chosen formula as an expression
print(f"The estimated weight change is given by the formula:\n{formula_str}\n")

# Calculate and print the numerical result by substituting values
# Breaking down the calculation for clarity
mass_flow_rate_str = f"({math.pi:.5f} * {d}**2 * {h} * {rho}) / (4 * {t})"
velocity_str = f"sqrt(2 * {g} * {H})"
full_calculation_str = f"Delta_W = {mass_flow_rate_str} * {velocity_str}"

mass_flow_rate_val = (math.pi * d**2 * h * rho) / (4 * t)
velocity_val = math.sqrt(2 * g * H)
delta_W = mass_flow_rate_val * velocity_val

print("Substituting the numerical values:")
print(full_calculation_str)
print(f"Delta_W = {mass_flow_rate_val:.4e} kg/s * {velocity_val:.4f} m/s")
print(f"Delta_W = {delta_W:.4e} N")
