import math

# --- Physical Constants ---
# Gravitational constant
G = 6.67430e-11  # m^3 kg^-1 s^-2
# Mass of the Sun
M_sun = 1.989e30  # kg
# Mass of the Earth
M_earth = 5.972e24  # kg
# Mass of the Moon
M_moon = 7.342e22  # kg
# Average distance from Earth to Sun (1 AU)
R_earth_sun = 1.496e11  # m
# Average distance from Moon to Earth
R_moon_earth = 3.844e8  # m
# Earth's orbital period
T_earth = 365.25 * 24 * 3600  # s
# Moon's orbital period
T_moon = 27.3 * 24 * 3600  # s

# --- 1. Calculate the ratio of Gravitational Forces (Sun's on Moon vs. Earth's on Moon) ---
# The distance from the Sun to the Moon is approximately the same as the Sun-Earth distance
# Force = G * M1 * M2 / r^2
F_sun_on_moon = (G * M_sun * M_moon) / (R_earth_sun**2)
F_earth_on_moon = (G * M_earth * M_moon) / (R_moon_earth**2)
force_ratio = F_sun_on_moon / F_earth_on_moon

print("--- Analysis of the Moon's Orbit ---")
print("\n1. Force Analysis:")
print(f"Force of Sun on Moon (N): {F_sun_on_moon:.3e}")
print(f"Force of Earth on Moon (N): {F_earth_on_moon:.3e}")
print(f"Ratio (Sun's Force / Earth's Force): {force_ratio:.2f}")
print("Conclusion 1: The Sun's gravitational pull on the Moon is more than twice as strong as the Earth's. Therefore, the Moon's trajectory is always curved towards the Sun. This rules out paths with inflection points (like D) or loops (like E).")


# --- 2. Calculate the ratio of Orbital Speeds ---
# Speed = 2 * pi * r / T
V_earth_around_sun = (2 * math.pi * R_earth_sun) / T_earth
V_moon_around_earth = (2 * math.pi * R_moon_earth) / T_moon
speed_ratio = V_earth_around_sun / V_moon_around_earth

print("\n2. Speed Analysis:")
print(f"Earth's speed around Sun (m/s): {V_earth_around_sun:,.0f}")
print(f"Moon's speed around Earth (m/s): {V_moon_around_earth:,.0f}")
print(f"Ratio (Earth's speed / Moon's speed): {speed_ratio:.2f}")
print("Conclusion 2: The Earth's speed orbiting the Sun is about 30 times faster than the Moon's speed orbiting the Earth. This means the Moon is always carried forward and never moves backward relative to the Sun. This confirms that the path cannot have loops (like E).")

print("\nFinal Conclusion:")
print("The Moon's orbit is essentially the Earth's large, curved path around the Sun (like B), but with a small 'wobble' due to its orbit around the Earth. Since the path is always convex towards the sun and has no loops, shape C is the correct representation.")