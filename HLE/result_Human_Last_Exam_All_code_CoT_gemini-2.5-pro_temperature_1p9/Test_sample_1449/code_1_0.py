import math

# Given parameters
m = 100 / 1000 # Mass of the block in kg
M = 10.0 # Mass of the wedge in kg
theta_deg = 30.0 # Angle of the wedge in degrees
h = 2.0 # Height of the wedge in m
mu = 1.0/2.0 # Coefficient of friction
g = 10.0 # Acceleration due to gravity in m/s^2

# Convert angle to radians for use in math functions
theta_rad = math.radians(theta_deg)
s = math.sin(theta_rad)
c = math.cos(theta_rad)

print("Step 1: Calculate the acceleration of the block relative to the wedge (a_rel).")
print("The formula for relative acceleration is a_rel = Num / Den, where:")
print("Num = g * (M + m) * (sin(theta) - mu * cos(theta))")
print("Den = M + m * sin(theta)^2 - m * mu * sin(theta) * cos(theta)")
print("-" * 20)

# Calculate the numerator and denominator of the acceleration formula
numerator_a = g * (M + m) * (s - mu * c)
denominator_a = M + m * s**2 - m * mu * s * c

print(f"Values used:")
print(f"g={g}, M={M}, m={m}, theta={theta_deg} deg, mu={mu}")
print(f"sin({theta_deg}) = {s:.4f}")
print(f"cos({theta_deg}) = {c:.4f}")
print("")

print("Calculation of the terms:")
print(f"Num = {g} * ({M} + {m}) * ({s:.4f} - {mu} * {c:.4f}) = {numerator_a:.4f}")
print(f"Den = {M} + {m} * {s:.4f}^2 - {m} * {mu} * {s:.4f} * {c:.4f} = {denominator_a:.4f}")
a_rel = numerator_a / denominator_a
print(f"a_rel = {numerator_a:.4f} / {denominator_a:.4f} = {a_rel:.4f} m/s^2")
print("\n")


print("Step 2: Calculate the distance the block slides along the wedge (d).")
print("The formula is d = h / sin(theta)")
print("-" * 20)
d = h / s
print(f"d = {h} / {s:.4f} = {d:.4f} m")
print("\n")

print("Step 3: Calculate the time (t) it takes for the block to slide down.")
print("Using the kinematic equation d = 0.5 * a_rel * t^2, we solve for t:")
print("t = sqrt(2 * d / a_rel)")
print("-" * 20)
t_squared = 2 * d / a_rel
t = math.sqrt(t_squared)
print(f"t = sqrt(2 * {d:.4f} / {a_rel:.4f})")
print(f"t = sqrt({t_squared:.4f})")
print(f"t = {t:.4f} s")
print("-" * 20)

print(f"\nThe exact amount of time it takes for the block to slide to the bottom is {t} seconds.")

<<<3.4392>>>