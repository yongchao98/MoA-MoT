import math

# Define the constants and given values from the problem statement.
# Gravitational constant in m^3 kg^-1 s^-2
G = 6.67 * 10**-11
# Solar mass in kg
M_solar = 1.99 * 10**30
# Side of the equilateral triangle in meters
L = 1.2 * 10**10
# Tangential velocity in km/s, converted to m/s
v = 125 * 1000

# The physical model equates the net gravitational force on a star
# with the centripetal force required for its orbit. This simplifies to:
# m = (v^2 * L) / G
# We will use this equation to find the mass of a single star.

# Calculate the mass of one star in kilograms
mass_in_kg = (v**2 * L) / G

# Convert the mass from kilograms to solar masses
mass_in_solar_masses = mass_in_kg / M_solar

# Print the final equation with the numbers plugged in, as requested.
print("The simplified equation for the mass (m) of a single component is:")
print("m = (v^2 * L) / G")
print("\nPlugging in the values:")
print(f"m = ({v} m/s)^2 * ({L} m) / ({G} m^3 kg^-1 s^-2)")

# Print the final answer rounded to one decimal place.
print(f"\nThe mass of a single component is {round(mass_in_solar_masses, 1)} solar masses.")

print(f"\n<<<{round(mass_in_solar_masses, 1)}>>>")