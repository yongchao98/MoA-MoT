import numpy as np

# Physical constants
M_sun = 1.989e30  # Mass of the Sun in kg
M_earth = 5.972e24 # Mass of the Earth in kg
R_sun_earth = 1.496e11 # Average distance from the Sun to the Earth in meters
r_earth_moon = 3.844e8  # Average distance from the Earth to the Moon in meters

# The distance from the Sun to the Moon is approximately R_sun_earth.
# The gravitational force is F = G * m1 * m2 / r^2.
# We want to find the ratio: F_sun_moon / F_earth_moon
# F_sun_moon = G * M_sun * M_moon / R_sun_earth^2
# F_earth_moon = G * M_earth * M_moon / r_earth_moon^2
# The ratio simplifies to (M_sun / M_earth) * (r_earth_moon / R_sun_earth)^2

print("Step 1: Define the physical constants.")
print(f"Mass of the Sun (M_sun): {M_sun:.3e} kg")
print(f"Mass of the Earth (M_earth): {M_earth:.3e} kg")
print(f"Sun-Earth distance (R_sun_earth): {R_sun_earth:.3e} m")
print(f"Earth-Moon distance (r_earth_moon): {r_earth_moon:.3e} m\n")

print("Step 2: Calculate the ratio of forces.")
print("Ratio = (M_sun / M_earth) * (r_earth_moon / R_sun_earth)^2\n")

# Calculate the components of the ratio
mass_ratio = M_sun / M_earth
distance_ratio_sq = (r_earth_moon / R_sun_earth)**2
total_ratio = mass_ratio * distance_ratio_sq

print("Step 3: Print the final result and interpret.")
print(f"The ratio of masses (M_sun / M_earth) is: {mass_ratio:.2f}")
print(f"The ratio of distances squared (r_earth_moon / R_sun_earth)^2 is: {distance_ratio_sq:.2e}")
print(f"The final ratio of forces (Sun's pull / Earth's pull) is:")
print(f"({M_sun:.3e} / {M_earth:.3e}) * ({r_earth_moon:.3e} / {R_sun_earth:.3e})^2 = {total_ratio:.2f}\n")

print("Interpretation:")
print(f"The Sun's gravitational pull on the Moon is approximately {total_ratio:.2f} times stronger than the Earth's pull.")
print("Because the Sun's force is always dominant, the Moon's orbit is always concave towards the Sun.")
print("The path never loops or curves away from the Sun. It is essentially the Earth's path around the Sun with a small monthly 'wobble' superimposed on it.")
print("This corresponds to option C, which shows a path that is always curved in one direction but has a periodic wave.")
