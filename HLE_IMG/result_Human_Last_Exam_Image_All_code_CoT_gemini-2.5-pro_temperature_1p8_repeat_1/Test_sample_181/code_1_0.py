import numpy as np

# Step 1: Define physical constants
G = 6.67430e-11  # Gravitational constant in m^3 kg^-1 s^-2
M_sun = 1.989e30  # Mass of the Sun in kg
M_earth = 5.972e24 # Mass of the Earth in kg
M_moon = 7.342e22  # Mass of the Moon in kg
R_sun_earth = 1.496e11 # Average distance from Sun to Earth in meters
r_earth_moon = 3.844e8   # Average distance from Earth to Moon in meters

# Step 2: Calculate the gravitational force of the Sun on the Moon
# We approximate the Sun-Moon distance as the Sun-Earth distance
F_sun_moon = G * (M_sun * M_moon) / (R_sun_earth**2)

# Step 3: Calculate the gravitational force of the Earth on the Moon
F_earth_moon = G * (M_earth * M_moon) / (r_earth_moon**2)

# Step 4: Compare the two forces
ratio = F_sun_moon / F_earth_moon

# Step 5: Print the results and explain the conclusion
print("Analysis of forces on the Moon:")
print("-" * 40)
print(f"Force of the Sun on the Moon (F_sun_moon): {F_sun_moon:.2e} N")
print(f"Force of the Earth on the Moon (F_earth_moon): {F_earth_moon:.2e} N")
print("-" * 40)
print(f"The force from the Sun is {ratio:.2f} times stronger than the force from the Earth.")
print("\nConclusion:")
print("Because the Sun's gravitational pull on the Moon is more than twice as strong as the Earth's pull, the Moon's trajectory is primarily an orbit around the Sun.")
print("The Earth's gravity only causes a small wobble in this orbit. Therefore, the Moon's path is always curved (concave) towards the Sun and never loops back on itself.")
print("This corresponds to diagram C, which shows a large curve with a smaller oscillation superimposed on it, without ever changing the direction of its main curvature.")
