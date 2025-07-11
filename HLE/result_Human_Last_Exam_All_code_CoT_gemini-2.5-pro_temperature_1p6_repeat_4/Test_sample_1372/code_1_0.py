import math

# --- Given Data ---
# Side of the equilateral triangle (m)
L = 1.2 * 10**10
# Tangential velocity (m/s). 1 km = 1000 m
v = 125 * 1000
# Gravitational constant (kg^-1 m^3 s^-2)
G = 6.67 * 10**-11
# Solar mass (kg)
M_solar = 1.99 * 10**30

# The equation to find the mass M of a single star is derived from
# equating the net gravitational force with the centripetal force.
# The simplified equation is M = (v^2 * L) / G.

print("Equation to calculate the mass M of a single component:")
print("M = (v^2 * L) / G\n")

print("Plugging in the given values:")
# We use scientific notation for G to maintain clarity with the problem statement.
print(f"M = ({v:.0f}^2 * {L:.1e}) / {G:.2e}")

# Calculate the mass of one star in kilograms
mass_kg = (v**2 * L) / G

print(f"\nResulting mass in kilograms: {mass_kg:.4e} kg")

# Convert the mass to solar masses
mass_solar_masses = mass_kg / M_solar

print(f"\nConverting to solar masses (1 solar mass = {M_solar:.2e} kg):")
print(f"Mass = {mass_kg:.4e} / {M_solar:.2e}")
print(f"\nThe mass of a single component is {mass_solar_masses:.1f} solar masses.")

<<<1.4>>>