import math

# Step 1: Define physical constants
G = 6.67430e-11  # Gravitational constant in m^3 kg^-1 s^-2
M_sun = 1.989e30      # Mass of the Sun in kg
M_earth = 5.972e24    # Mass of the Earth in kg
M_moon = 7.342e22     # Mass of the Moon in kg

# Average distances in meters
R_earth_sun = 1.496e11  # Earth-Sun distance (1 AU)
R_earth_moon = 3.844e8   # Earth-Moon distance

# The distance from the Sun to the Moon varies, but we can approximate it
# with the Earth-Sun distance for this magnitude comparison.
R_sun_moon_approx = R_earth_sun

# Step 2: Calculate gravitational forces
# Force of the Sun on the Moon
F_sun_moon = (G * M_sun * M_moon) / (R_sun_moon_approx**2)

# Force of the Earth on the Moon
F_earth_moon = (G * M_earth * M_moon) / (R_earth_moon**2)

# Step 3: Compare the forces
ratio = F_sun_moon / F_earth_moon

print(f"Force of the Sun on the Moon: {F_sun_moon:.2e} N")
print(f"Force of the Earth on the Moon: {F_earth_moon:.2e} N")
print(f"Ratio (Sun's force / Earth's force): {ratio:.2f}")

print("\nConclusion:")
print("The Sun's gravitational force on the Moon is more than twice as strong as the Earth's.")
print("This means the Moon's trajectory is primarily an orbit around the Sun, which is perturbed by the Earth.")
print("The path is therefore always concave towards the Sun, as shown in diagram C. It wobbles but does not create loops.")
