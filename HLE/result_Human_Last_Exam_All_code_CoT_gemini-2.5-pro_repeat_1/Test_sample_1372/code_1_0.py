import math

# --- Given Data ---
# Gravitational constant in kg^-1 m^3 s^-2
G = 6.67 * 10**-11
# Solar mass in kg
M_sun = 1.99 * 10**30
# Side length of the equilateral triangle in m
L = 1.2 * 10**10
# Tangential velocity in km/s
v_kms = 125

# --- Conversion ---
# Convert velocity to m/s
v_ms = v_kms * 1000

# --- Derivation Summary ---
# The net gravitational force on one star is F_grav = sqrt(3) * G * m^2 / L^2.
# The centripetal force is F_cent = m * v^2 / r, where the orbital radius r = L / sqrt(3).
# Equating F_grav and F_cent and simplifying gives the formula for mass m.

# --- Calculation ---
# m = v^2 * L / G
# Mass of a single star in kg
m_kg = (v_ms**2 * L) / G

# Mass of a single star in solar masses
m_solar = m_kg / M_sun

# --- Output ---
print("Calculating the mass of a single star component.")
print("The mass (m) can be found by equating the net gravitational force with the centripetal force.")
print("The resulting equation is: m = v^2 * L / G")
print("\nUsing the following values:")
print(f"v (velocity) = {v_kms} km/s = {v_ms:.2e} m/s")
print(f"L (side length) = {L:.2e} m")
print(f"G (gravitational constant) = {G:.2e} kg^-1 m^3 s^-2")
print("\nPlugging the numbers into the equation:")
print(f"m = ({v_ms:.2e})^2 * ({L:.2e}) / ({G:.2e})")

# Calculate intermediate steps for clarity
v_squared = v_ms**2
v_squared_L = v_squared * L
print(f"m = ({v_squared:.4e}) * ({L:.2e}) / ({G:.2e})")
print(f"m = ({v_squared_L:.4e}) / ({G:.2e})")
print(f"m = {m_kg:.4e} kg")

print("\nTo convert the mass from kg to solar masses, divide by the solar mass:")
print(f"Mass in solar masses = {m_kg:.4e} kg / {M_sun:.2e} kg")
print(f"Mass of a single component = {m_solar:.1f} solar masses")
<<<1.4>>>