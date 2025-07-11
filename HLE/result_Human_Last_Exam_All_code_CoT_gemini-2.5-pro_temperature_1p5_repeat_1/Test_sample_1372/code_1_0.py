import math

# --- Given Data and Constants ---
# Gravitational constant in m^3 kg^-1 s^-2
G = 6.67 * 10**-11
# Solar mass in kg
solar_mass = 1.99 * 10**30
# Side of the equilateral triangle in m
a = 1.2 * 10**10
# Tangential velocity in km/s
v_kms = 125

# --- Calculations ---
# Convert velocity from km/s to m/s
v_ms = v_kms * 1000

# Based on the force balance, the mass M of a single star can be calculated using the formula: M = v^2 * a / G
# Calculate the mass of a single component in kg
mass_kg = (v_ms**2 * a) / G

# Convert the mass from kg to solar masses
mass_solar = mass_kg / solar_mass

# --- Output ---
print("The formula for the mass (M) of one star is M = v^2 * a / G")
print("Substituting the values into the equation:")
print(f"M = ({v_ms:.0f} m/s)^2 * ({a:.1e} m) / ({G:.2e} m^3 kg^-1 s^-2)")
print(f"Mass in kg = {mass_kg:.4e} kg")
print("\nConverting mass to solar masses:")
print(f"Mass (solar masses) = {mass_kg:.4e} kg / {solar_mass:.2e} kg/solar_mass")
print(f"\nThe mass of a single component is {mass_solar:.1f} solar masses.")
