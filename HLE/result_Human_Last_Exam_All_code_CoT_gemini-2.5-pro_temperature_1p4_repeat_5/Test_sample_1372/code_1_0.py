import math

# --- Data and Constants ---

# Side of the equilateral triangle (m)
s = 1.2 * 10**10

# Tangential velocity (m/s). Given as 125 km/s.
v = 125 * 1000

# Gravitational constant (kg^-1 m^3 s^-2)
G = 6.67 * 10**-11

# Solar mass (kg)
M_solar = 1.99 * 10**30

# The provided period of 1.9 days is inconsistent with the given tangential velocity
# for this stable configuration. The calculation will proceed using the velocity.

# --- Calculation ---

# The formula for the mass of a single component is derived by equating
# the net gravitational force with the centripetal force: m = v^2 * s / G.
# First, calculate the mass in kilograms.
m_kg = (v**2 * s) / G

# Second, convert the mass to solar masses.
m_solar_masses = m_kg / M_solar

# --- Output ---

# As requested, output the final equation with the numerical values.
# The equation calculates the mass in solar masses directly.
print("Final Calculation:")
print(f"Mass_in_solar_masses = ( (Velocity)^2 * Side_Length / Gravitational_Constant ) / Solar_Mass")
print(f"Mass_in_solar_masses = ( ({v:.2e})^2 * {s:.2e} / {G:.2e} ) / {M_solar:.2e}")
print("")
# Print the final result, formatted to one decimal place.
print(f"The mass of a single component is {m_solar_masses:.1f} solar masses.")
