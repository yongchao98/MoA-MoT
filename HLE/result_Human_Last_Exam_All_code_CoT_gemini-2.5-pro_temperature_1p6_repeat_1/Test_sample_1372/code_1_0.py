import math

# Define the given constants and data
G = 6.67 * 10**-11  # Gravitational constant in kg^-1 m^3 s^-2
L = 1.2 * 10**10     # Side of the equilateral triangle in meters
v_kms = 125          # Tangential velocity in km/s
M_sun = 1.99 * 10**30 # Solar mass in kg

# Convert velocity from km/s to m/s
v_ms = v_kms * 1000

# Calculate the mass of a single star in kg using the derived formula: m = v^2 * L / G
mass_kg = (v_ms**2 * L) / G

# Convert the mass to solar masses
mass_solar = mass_kg / M_sun

# Print the final equation with numerical values
print(f"The equation for the mass (m) of a single component is:")
print(f"m = (v^2 * L) / G")
print(f"Substituting the values:")
print(f"m = ({v_ms:.0f} m/s)^2 * ({L:.1e} m) / ({G:.2e} kg^-1 m^3 s^-2)")

# Print the result in solar masses
print(f"\nThe mass of a single component star is {mass_solar:.1f} solar masses.")

# Final answer in the specified format
# The formatted output string is sufficient, but let's be explicit
# print(f'<<<{mass_solar:.1f}>>>')