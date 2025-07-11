import math

# --- Step 1: Define Constants ---
h = 1.0  # robot height in meters
r = 0.25  # arm length in meters
l_c = 10.0  # chain length in meters
v = 10.0  # robot speed in m/s
dot_beta = 1.0  # arm angular speed in rad/s
d = 20.0  # visible path length in meters
l_shadow = 10 * math.sqrt(3)  # shadow length in meters

# --- Step 2: Calculate Derived Parameters ---
# Radius of the circular path
R = d / 2.0
# Angular speed of the robot on the path
omega = v / R  # This is 1.0 rad/s, same as dot_beta

# --- Step 3: Set up the Quadratic Equation for sin(t) ---
# The condition for the chain tip to be at height l_c leads to a
# quadratic equation of the form a*x^2 + b*x + c = 0, where x = sin(t).
# Based on the derivation in the thinking steps:
# a*sin(t)^2 + b*sin(t) + c = 0
# with sin(t)^2 - (40 + sqrt(3))*sin(t) + (39 - 4*sqrt(3)) = 0

a = 1.0
b = -(40.0 + math.sqrt(3))
c = 39.0 - 4.0 * math.sqrt(3)

# --- Step 4: Output the equation ---
print("The problem reduces to solving a quadratic equation for sin(t):")
# Use '+' for negative b and c to format the equation string correctly
b_sign = "-" if b < 0 else "+"
c_sign = "-" if c < 0 else "+"
print(f"{a:.3f}*sin(t)^2 {b_sign} {abs(b):.3f}*sin(t) {c_sign} {abs(c):.3f} = 0")
print("-" * 20)


# --- Step 5: Solve the Quadratic Equation ---
# Calculate the discriminant
discriminant = b**2 - 4 * a * c

if discriminant < 0:
    print("No real solution for t exists.")
else:
    # Calculate the two possible solutions for sin(t)
    sin_t1 = (-b + math.sqrt(discriminant)) / (2 * a)
    sin_t2 = (-b - math.sqrt(discriminant)) / (2 * a)

    # The value of sin(t) must be between -1 and 1.
    # We choose the valid root.
    valid_sin_t = None
    if -1 <= sin_t1 <= 1:
        valid_sin_t = sin_t1
    elif -1 <= sin_t2 <= 1:
        valid_sin_t = sin_t2

    if valid_sin_t is not None:
        # Calculate the time t. We want the first positive time, so we take the principal value of arcsin.
        t = math.asin(valid_sin_t)
        print(f"The first time the chain loses contact with the ground is at t = {t:.4f} seconds.")
        print(f"<<<{t:.4f}>>>")
    else:
        print("No valid solution for sin(t) found in the range [-1, 1].")
