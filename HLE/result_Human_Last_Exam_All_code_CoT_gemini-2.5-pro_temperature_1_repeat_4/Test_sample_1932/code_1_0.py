import math

# Given parameters
d = 0.01  # diameter in meters
h = 0.02  # height of sand column in meters
H = 0.04  # height of the chamber in meters
rho = 1500  # density in kg/m^3
t = 60  # time in seconds
g = 9.8  # acceleration due to gravity in m/s^2

# Expression for the change in weight (Option B)
# ΔW = - (π * d² * h * ρ / (4 * t)) * sqrt(2 * g * H)

# Calculate the mass flow rate (μ)
mu_numerator = math.pi * d**2 * h * rho
mu_denominator = 4 * t
mu = mu_numerator / mu_denominator

# Calculate the impact velocity at the start (v)
v = math.sqrt(2 * g * H)

# Calculate the change in weight
delta_W = -mu * v

# Print the calculation steps with numbers
print("Step 1: Calculate the mass flow rate (μ)")
print(f"μ = (π * d² * h * ρ) / (4 * t)")
print(f"μ = (π * ({d})² * {h} * {rho}) / (4 * {t})")
print(f"μ = {mu:.4e} kg/s\n")

print("Step 2: Calculate the fall velocity at the start (v)")
print(f"v = sqrt(2 * g * H)")
print(f"v = sqrt(2 * {g} * {H})")
print(f"v = {v:.4f} m/s\n")

print("Step 3: Calculate the change in weight (ΔW)")
print("ΔW = -μ * v")
print(f"ΔW = -({mu:.4e}) * {v:.4f}")
print(f"ΔW = {delta_W:.4e} N")
