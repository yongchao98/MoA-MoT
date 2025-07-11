import math

# Given data
G = 6.67 * 10**-11  # Gravitational constant in kg^-1 m^3 s^-2
SOLAR_MASS = 1.99 * 10**30  # Solar mass in kg
a = 1.2 * 10**10  # Side of the equilateral triangle in meters
v_km_s = 125  # Tangential velocity in km/s

# Convert velocity to m/s
v_m_s = v_km_s * 1000

# The equation for the mass of a single component is derived from F_gravity = F_centripetal
# (G * m^2 / a^2) * sqrt(3) = m * v^2 / (a / sqrt(3))
# G * m / a = v^2
# m = a * v^2 / G

# Calculate the mass of a single star in kg
mass_kg = (a * v_m_s**2) / G

# Convert the mass to solar masses
mass_solar = mass_kg / SOLAR_MASS

# Print the final equation with values
print(f"The equation to solve for the mass 'm' in kg is: m = (a * v^2) / G")
print(f"m = ({a} * {v_m_s}^2) / {G}")
print(f"Mass of one star: {mass_kg:.3e} kg")

# Print the final answer in solar masses, rounded to one decimal place
print(f"Mass in solar masses: {mass_solar:.1f}")

# Final answer in the required format
final_answer = round(mass_solar, 1)
# print(f'<<<{final_answer}>>>') # This is for internal check, not final output.
# The prompt wants the answer at the end, outside the code block.