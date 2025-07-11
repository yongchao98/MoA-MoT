import math

# Description:
# This script calculates the mass of a single star in a symmetric triple star system.
# The calculation is based on the principle that the gravitational force on a star
# provides the necessary centripetal force for its orbit.

# Given data from the problem
side_length_L = 1.2e10  # Side of the equilateral triangle in meters
velocity_v_kms = 125    # Tangential velocity in km/s

# Physical constants
GRAVITATIONAL_CONSTANT_G = 6.67e-11 # m^3 kg^-1 s^-2
SOLAR_MASS = 1.99e30                 # kg

# Convert velocity from km/s to m/s for SI unit consistency
velocity_v_ms = velocity_v_kms * 1000

# The derived formula for the mass 'm' of a single component is m = (v^2 * L) / G.
# Let's print the equation with the numbers substituted in.
print("Derived equation for mass (m):")
print("m = v^2 * L / G")
print("\nSubstituting the given values:")
print(f"m = ({velocity_v_ms:.3e} m/s)^2 * ({side_length_L:.3e} m) / ({GRAVITATIONAL_CONSTANT_G:.3e} m^3 kg^-1 s^-2)")

# Calculate the mass of a single star in kilograms
mass_kg = (velocity_v_ms**2 * side_length_L) / GRAVITATIONAL_CONSTANT_G

# Convert the mass from kilograms to solar masses
mass_solar_masses = mass_kg / SOLAR_MASS

# Print the final answer rounded to one decimal place
print(f"\nThe mass of a single component is {mass_solar_masses:.1f} solar masses.")
