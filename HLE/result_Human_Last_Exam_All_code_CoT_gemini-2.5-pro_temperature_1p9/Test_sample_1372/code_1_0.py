import math

# Define the given constants and data
G = 6.67e-11      # Gravitational constant in m^3 kg^-1 s^-2
L = 1.2e10        # Side of the equilateral triangle in meters
v_kms = 125       # Tangential velocity in km/s
solar_mass = 1.99e30  # Mass of the sun in kg

# Convert velocity from km/s to m/s
v_ms = v_kms * 1000

# Calculate the mass of a single star in kg using the formula m = v^2 * L / G
mass_kg = (v_ms**2 * L) / G

# Convert the mass to solar masses
mass_solar = mass_kg / solar_mass

# Print the final equation with all values substituted
print("The mass 'm' is calculated using the formula: m = v^2 * L / G")
print("Substituting the values into the equation:")
print(f"m = ({v_ms:.0f} m/s)^2 * ({L:.1e} m) / ({G:.2e} m^3 kg^-1 s^-2)")
print(f"Calculated mass in kg: {mass_kg:.3e} kg")
print(f"Calculated mass in solar masses: {mass_solar:.3f} M_sun")
print(f"The mass of a single component, accurate to one decimal place, is {mass_solar:.1f} solar masses.")
