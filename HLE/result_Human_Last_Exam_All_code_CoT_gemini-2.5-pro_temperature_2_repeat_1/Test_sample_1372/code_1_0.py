import math

# Constants
G = 6.67 * 10**-11  # Gravitational constant in kg^-1 m^3 s^-2
M_sun = 1.99 * 10**30 # Solar mass in kg

# Given data
a = 1.2 * 10**10     # Side of the equilateral triangle in meters
v_kms = 125          # Tangential velocity in km/s
v_ms = v_kms * 1000  # Convert velocity to m/s

# The derived formula to calculate the mass of a single component is:
# m = (v^2 * a) / G

# Calculate the mass in kg
mass_kg = (v_ms**2 * a) / G

# Convert the mass to solar masses
mass_solar = mass_kg / M_sun

# Print the final equation with all the numbers
print("Equation to calculate mass in kg:")
print(f"m = ({v_ms:.2e}^2 * {a:.2e}) / {G:.2e}")

# Print the intermediate and final results
print(f"\nMass of a single component: {mass_kg:.4e} kg")
print(f"Mass in solar masses: {mass_solar:.1f} M_sun")

# Final answer in the required format
final_answer = round(mass_solar, 1)
# The line below is for the final submission format.
# print(f"<<<{final_answer}>>>")