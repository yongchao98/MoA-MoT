import math

# Step 1: Define the given parameters
h = 1.0  # m, robot height
r = 0.25  # m, arm length
l_c = 10.0  # m, chain length
v = 10.0  # m/s, robot speed
d = 20.0  # m, visible diameter of the path
l_shadow = 10.0 * math.sqrt(3)  # m, shadow length

# Step 2: Calculate the geometry of the path
R = d / 2.0  # Radius of the circular path
# The tilt angle alpha is found from the shadow projection
# l_shadow = d * cos(alpha)
cos_alpha = l_shadow / d
alpha = math.acos(cos_alpha)
sin_alpha = math.sin(alpha)

# Step 3: Formulate the equation for the height of the arm's tip
# The height of the arm tip z_tip(t) must equal the chain length l_c
# for the chain to lose contact with the ground.
# z_tip(t) = z_feet(t) + z_body_vertical + z_arm_vertical
# The robot starts at theta_0 = pi/2 and moves with angular speed omega=v/R=1 rad/s.
# So, the angle at time t is theta(t) = pi/2 + t.
# z_feet(t) = R * sin(alpha) * (1 - cos(theta(t))) = R * sin(alpha) * (1 + sin(t))
# z_body_vertical = h * cos(alpha)
# z_arm_vertical(t) = r * (sin(alpha) * cos^2(t) + cos(alpha) * sin(t))
# The equation is:
# R*sin(alpha)*(1+sin(t)) + h*cos(alpha) + r*(sin(alpha)*cos^2(t) + cos(alpha)*sin(t)) = l_c
# Let s = sin(t), then cos^2(t) = 1 - s^2.
# R*sin(a)*(1+s) + h*cos(a) + r*(sin(a)*(1-s^2) + cos(a)*s) - l_c = 0
# Rearranging into a quadratic equation for s: As^2 + Bs + C = 0
# A = -r * sin(alpha)
# B = R * sin(alpha) + r * cos(alpha)
# C = R * sin(alpha) + h * cos(alpha) + r * sin(alpha) - l_c

A = -r * sin_alpha
B = R * sin_alpha + r * cos_alpha
C = R * sin_alpha + h * cos_alpha + r * sin_alpha - l_c

print("The problem reduces to solving for t in the equation z_tip(t) = l_c.")
print("The equation, with all values substituted, is:")
# Print each term of the equation before grouping them into a quadratic
# Equation: term1 * (1 + sin(t)) + term2 + term3 * (1-sin(t)^2) + term4 * sin(t) - l_c = 0
term1 = R * sin_alpha
term2 = h * cos_alpha
term3 = r * sin_alpha
term4 = r * cos_alpha
print(f"{term1:.4f} * (1 + sin(t)) + {term2:.4f} + {term3:.4f} * cos^2(t) + {term4:.4f} * sin(t) - {l_c:.4f} = 0")
print("\nThis simplifies to the quadratic equation A*s^2 + B*s + C = 0, where s = sin(t):")
print(f"{A:.4f}*s^2 + {B:.4f}*s + {C:.4f} = 0")

# Step 4: Solve the quadratic equation for s = sin(t)
discriminant = B**2 - 4 * A * C

# There are two potential solutions for s
s1 = (-B + math.sqrt(discriminant)) / (2 * A)
s2 = (-B - math.sqrt(discriminant)) / (2 * A)

# Choose the physically valid solution for s (-1 <= s <= 1)
if -1 <= s1 <= 1:
    s_sol = s1
else:
    s_sol = s2

# Step 5: Calculate the time t
# We need the first positive time, so t = arcsin(s)
t = math.asin(s_sol)

print(f"\nThe solution for sin(t) is {s_sol:.4f}.")
print(f"The first time t when the chain loses contact with the ground is:")
print(f"{t:.4f} seconds")

# Final answer in the required format
print(f"\n<<<{t:.4f}>>>")