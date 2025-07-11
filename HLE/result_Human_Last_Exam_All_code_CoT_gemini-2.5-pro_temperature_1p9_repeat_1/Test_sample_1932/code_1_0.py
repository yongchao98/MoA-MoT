import math

# Given parameters
d_cm = 1  # diameter in cm
h_cm = 2  # height of sand column in cm
H_cm = 4  # height of one chamber in cm
rho_kg_m3 = 1500  # density in kg/m^3
t_min = 1  # time in minutes

# Constants
g = 9.8  # m/s^2

# Convert units to SI
d = d_cm / 100.0  # m
h = h_cm / 100.0  # m
H = H_cm / 100.0  # m
t = t_min * 60.0   # s

# The chosen expression for the weight change from option B
# Delta_W = - (mass_flow_rate) * (fall_velocity)
# Delta_W = - ( (pi * d**2 * h * rho) / (4 * t) ) * sqrt(2 * g * H)

# Calculate each part of the formula
mass_flow_rate = (math.pi * d**2 * h * rho_kg_m3) / (4 * t)
fall_velocity = math.sqrt(2 * g * H)
delta_W = -mass_flow_rate * fall_velocity

# Print the explanation and the formula with values
print("The analysis suggests the hourglass becomes lighter while running.")
print("The estimated change in weight, Delta W, is given by the formula from option B:")
print("Delta_W = - (pi * d^2 * h * rho / (4 * t)) * sqrt(2 * g * H)\n")
print("Plugging in the given values:")
print(f"d = {d:.2f} m")
print(f"h = {h:.2f} m")
print(f"H = {H:.2f} m")
print(f"rho = {rho_kg_m3} kg/m^3")
print(f"t = {t} s")
print(f"g = {g} m/s^2\n")

# Output the equation with numbers plugged in, as requested.
# Using unicode characters for better readability
print("The final equation is:")
print(f"ΔW = - (π * ({d:.2f})² * {h:.2f} * {rho_kg_m3} / (4 * {t:.0f})) * √(2 * {g} * {H:.2f})")
print(f"   = - ({mass_flow_rate:.3e} kg/s) * ({fall_velocity:.3f} m/s)")
print(f"   = {delta_W:.3e} N")
