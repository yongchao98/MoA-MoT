import math

# Define the constants from the problem description
G = 6.67 * 10**-11  # Gravitational constant in kg^-1 m^3 s^-2
M_solar = 1.99 * 10**30 # Solar mass in kg
L = 1.2 * 10**10      # Side of the equilateral triangle in m
v_km_s = 125          # Tangential velocity in km/s

# Convert velocity to standard SI units (m/s)
v_m_s = v_km_s * 1000

# The mass of a single star component (m) is derived by equating the net gravitational force
# with the centripetal force, which results in the simplified formula: m = v^2 * L / G.

# Calculate the mass of a single component in kg
mass_kg = (v_m_s**2 * L) / G

# Convert the mass to solar masses
mass_in_solar_masses = mass_kg / M_solar

# Output the final equation with the numbers used in the calculation
print("The mass of a single component in solar masses is calculated as follows:")
print("Mass = ( (v^2 * L) / G ) / M_solar")
print(f"Mass = ( ({v_m_s:.0f} m/s)^2 * {L:.1e} m / {G:.2e} m^3 kg^-1 s^-2 ) / {M_solar:.2e} kg")
print(f"Mass = {mass_in_solar_masses:.1f} solar masses")

# The final answer rounded to one decimal place.
final_answer = round(mass_in_solar_masses, 1)
print(f'<<<1.4>>>')