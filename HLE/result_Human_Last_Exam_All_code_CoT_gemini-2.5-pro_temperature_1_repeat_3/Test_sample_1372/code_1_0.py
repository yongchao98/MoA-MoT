import math

# Define the constants and given values
G = 6.67 * 10**-11  # Gravitational constant in m^3 kg^-1 s^-2
L = 1.2 * 10**10     # Side length of the triangle in meters
v_km_s = 125         # Tangential velocity in km/s
M_solar = 1.99 * 10**30 # Solar mass in kg

# Convert velocity to m/s
v = v_km_s * 1000

# The derived formula for the mass of a single component is m = (v^2 * L) / G
# Calculate the mass in kg
mass_kg = (v**2 * L) / G

# Convert the mass to solar masses
mass_solar = mass_kg / M_solar

# Print the equation with the numerical values
print("Equation for the mass of a single star:")
print(f"m = (v^2 * L) / G")
print(f"m = (({v:.2e} m/s)^2 * {L:.2e} m) / {G:.2e} m^3 kg^-1 s^-2")

# Print the final result in solar masses, rounded to one decimal place
print(f"\nThe mass of a single component is {mass_solar:.1f} solar masses.")

# Final answer in the required format
print(f"\n<<<{mass_solar:.1f}>>>")