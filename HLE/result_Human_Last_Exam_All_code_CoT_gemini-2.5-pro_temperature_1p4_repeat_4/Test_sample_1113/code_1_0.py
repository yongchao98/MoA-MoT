import math

# Given parameters
h = 1.0       # Robot height in m
r_arm = 0.25  # Arm length in m (25 cm)
l_c = 10.0    # Chain length in m
v = 10.0      # Robot speed in m/s
d = 20.0      # Visible diameter of the path in m
l_shadow = 10 * math.sqrt(3) # Shadow length in m

# Step 1: Determine path geometry
# The side-view length is the diameter of the circular path.
R = d / 2.0
# The tilt angle alpha is found from the diameter and shadow length (sun is overhead)
# l_shadow = d * cos(alpha)
cos_alpha = l_shadow / d
alpha = math.acos(cos_alpha)
sin_alpha = math.sin(alpha)

# Step 2: The height of the robot's base
# The robot reaches the highest point (z = d*sin(alpha) = 2*R*sin(alpha)) at t = (pi*R/2)/v = pi/2 s.
# The height of the base z_base(t) must satisfy this.
# A function z_base(t) = R*sin(alpha)*(1 + sin(t)) works.
# At t=pi/2, z_base = R*sin(alpha)*(1+1) = 2*R*sin(alpha).
# At t=3*pi/2 (lowest point), z_base = R*sin(alpha)*(1-1) = 0. This matches the problem statement.
# z_base(t) = R * sin(alpha) * (1 + sin(t))

# Step 3: Height of the arm's tip, z_tip(t)
# z_tip(t) is the sum of:
# 1. Robot base height: z_base(t)
# 2. Vertical component of robot's body: h_vert = h * cos(alpha)
# 3. Vertical component of arm's length: z_arm(t)
# The arm rotates with angle beta(t) = pi/2 + t in a plane spanned by the path's normal and tangent vectors.
# z_arm(t) = r_arm * (cos(t)**2 * sin_alpha - sin(t) * cos_alpha)
# z_tip(t) = z_base(t) + h_vert + z_arm(t)

# Step 4 & 5: Set up and solve the equation z_tip(t) = l_c
# The equation becomes a quadratic in x = sin(t): A*x^2 + B*x + C = 0
# (-r_arm*sin_alpha)*x^2 + (R*sin_alpha - r_arm*cos_alpha)*x + (R*sin_alpha + h*cos_alpha + r_arm*sin_alpha - l_c) = 0
A = -r_arm * sin_alpha
B = R * sin_alpha - r_arm * cos_alpha
C = R * sin_alpha + h * cos_alpha + r_arm * sin_alpha - l_c

# Solve the quadratic equation for x = sin(t)
discriminant = B**2 - 4 * A * C
# We expect two solutions for x
x1 = (-B + math.sqrt(discriminant)) / (2 * A)
x2 = (-B - math.sqrt(discriminant)) / (2 * A)

# The valid solution for sin(t) must be between -1 and 1
sin_t = 0
if -1 <= x1 <= 1:
    sin_t = x1
elif -1 <= x2 <= 1:
    sin_t = x2

# Find the first time t > 0
t = math.asin(sin_t)

# Print the explanation and the equation
print("The problem is solved by finding the time 't' when the height of the robot arm's tip equals the length of the chain.")
print("The equation for the height of the arm's tip, z_tip(t), is derived from the geometry of the system.")
print(f"Path radius R = {R:.1f} m, Path tilt angle alpha = {math.degrees(alpha):.1f} degrees.")
print("\nLet x = sin(t). The condition z_tip(t) = l_c can be written as a quadratic equation A*x^2 + B*x + C = 0.")
print("The coefficients are calculated as follows:")
print(f"A = -r_arm*sin(alpha) = -{r_arm}*{sin_alpha:.4f} = {A:.4f}")
print(f"B = R*sin(alpha) - r_arm*cos(alpha) = {R}*{sin_alpha:.4f} - {r_arm}*{cos_alpha:.4f} = {B:.4f}")
print(f"C = R*sin(alpha) + h*cos(alpha) + r_arm*sin(alpha) - l_c = {R}*{sin_alpha:.4f} + {h}*{cos_alpha:.4f} + {r_arm}*{sin_alpha:.4f} - {l_c} = {C:.4f}")

print("\nThe final equation to solve for sin(t) is:")
print(f"({A:.4f}) * sin(t)^2 + ({B:.4f}) * sin(t) + ({C:.4f}) = 0")
print(f"\nSolving this equation gives sin(t) = {sin_t:.4f}.")
print(f"The first time t > 0 when this occurs is t = asin({sin_t:.4f}).")
print(f"\nTime when the chain first loses contact with the ground: {t:.4f} seconds.")

# Final answer in the specified format
print(f"\n<<<{t:.4f}>>>")