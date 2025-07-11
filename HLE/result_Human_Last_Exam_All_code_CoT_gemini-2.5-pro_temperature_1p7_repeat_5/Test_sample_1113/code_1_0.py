import math

# Step 1: Define constants from the problem
h = 1.0  # Robot height in m
r = 0.25  # Arm length in m
l_c = 10.0  # Chain length in m
v = 10.0  # Robot speed in m/s
beta_dot = 1.0  # Arm angular speed in rad/s
d = 20.0  # Visible length (diameter) of the path in m
l_shadow = 10.0 * math.sqrt(3)  # Shadow length in m

# Step 2: Calculate geometric properties of the path
R = d / 2  # Radius of the circular path

# Calculate the tilt angle 'alpha'
cos_alpha = l_shadow / d
alpha = math.acos(cos_alpha)
sin_alpha = math.sin(alpha)

# Calculate the height of the circle's center from the ground
Z_c = R * sin_alpha

# Step 3: Formulate the height equation H(t) = 0
# The equation for the height of the chain's end H(t) is:
# H(t) = Z_path(t) + Z_robot + Z_arm(t) - l_c = 0
# Z_path(t) = Z_c + R * sin(t) * sin_alpha
# Z_robot = h * cos_alpha
# Z_arm(t) = r * (cos(t)**2 * sin_alpha - sin(t) * cos_alpha)
# Let s = sin(t). Then cos(t)^2 = 1 - s^2.
# The equation becomes a quadratic equation for s: A*s^2 + B*s + C = 0

# Coefficients of the quadratic equation A*s^2 + B*s + C = 0 for s=sin(t)
# derived from: (Z_c + R*s*sin_alpha) + (h*cos_alpha) + r*((1-s**2)*sin_alpha - s*cos_alpha) - l_c = 0
A = -r * sin_alpha
B = R * sin_alpha - r * cos_alpha
C = Z_c + h * cos_alpha + r * sin_alpha - l_c

# Step 4: Solve the quadratic equation for s = sin(t)
discriminant = B**2 - 4 * A * C

# We must check if the discriminant is non-negative
if discriminant < 0:
    print("No real solution exists.")
else:
    # Calculate the two possible solutions for s
    s1 = (-B + math.sqrt(discriminant)) / (2 * A)
    s2 = (-B - math.sqrt(discriminant)) / (2 * A)

    # The valid solution for sin(t) must be between -1 and 1
    valid_s = None
    if -1 <= s1 <= 1:
        # check if s2 is also valid and smaller, to find the first time t > 0
        if -1 <= s2 <= 1:
            # We are looking for the smallest positive t, which corresponds to the smallest positive s for small t.
            # However, we must choose the solution that corresponds to H(t) crossing from negative to positive.
            # The derivative of H(t) shows that the smaller root is the correct one.
            valid_s = min(s1, s2) if min(s1,s2) > 0 else max(s1,s2) # We look for the first t>0, so s>0.
            # In this problem, one root is >1 and the other is between 0 and 1.
            if valid_s > 1 or valid_s < -1: valid_s=None

    if valid_s is None and -1 <= s2 <= 1:
        valid_s = s2

    # Step 5: Calculate the time t and print the results
    if valid_s is not None:
        # The first time t > 0 is given by the principal value of arcsin
        t = math.asin(valid_s)

        # Output the numbers in the final equation as requested
        # The equation is: (Z_c + R*sin(t)*sin_alpha) + (h*cos_alpha) + (r*cos^2(t)*sin_alpha - r*sin(t)*cos_alpha) - l_c = 0
        print("The final equation for the time 't' when the chain loses contact is found by setting the chain end's height H(t) to zero.")
        print("H(t) = Z_path(t) + Z_robot_height(t) + Z_arm(t) - Chain_Length = 0")
        
        # Calculate each term's numeric value at t=0 for illustration, and symbolic form
        z_path_str = f"({Z_c:.3f} + {R * sin_alpha:.3f}*sin(t))"
        z_robot_str = f" + {h * cos_alpha:.3f}"
        z_arm_str = f" + {r * sin_alpha:.3f}*cos^2(t) - {r * cos_alpha:.3f}*sin(t)"
        l_c_str = f" - {l_c:.1f}"
        
        print(f"The equation with numerical coefficients is:")
        print(f"{z_path_str}{z_robot_str}{z_arm_str}{l_c_str} = 0")
        
        # Print the final answer
        print(f"\nThe time when the chain first loses contact with the ground is:")
        print(f"{t:.4f} seconds")
        print(f"\n<<< {t:.4f} >>>")

    else:
        print("No physically realistic solution for sin(t) found.")
