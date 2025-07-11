import math

# PLAN:
# 1. Define all given physical constants from the problem description.
# 2. Calculate the geometric properties of the tilted circular path: its radius R,
#    tilt angle alpha, and the height of its center from the ground.
# 3. Determine the robot's starting position on the circle based on the timing constraint
#    for reaching the highest point.
# 4. Formulate an equation for the vertical height of the tip of the robot's arm.
#    The chain loses contact when this height equals the chain's length.
# 5. Rearrange this equation into a quadratic form A*s^2 + B*s + C = 0, where s = sin(t).
# 6. Solve the quadratic equation for 's' and then find the smallest positive time 't'
#    from the valid solution(s).
# 7. Print the coefficients of the equation and the final calculated time.

# Step 1: Define constants
h = 1       # robot height in meters
r = 0.25    # arm length in meters
l_c = 10    # chain length in meters
v = 10      # robot speed in m/s
d = 20      # visible path length (diameter) in meters
l_shadow = 10 * math.sqrt(3) # shadow length in meters
dot_beta = 1 # arm angular speed in rad/s

# Step 2: Calculate geometric parameters
R = d / 2
cos_alpha = l_shadow / d
alpha = math.acos(cos_alpha)
sin_alpha = math.sin(alpha)
h_center = R * sin_alpha

# Step 3 & 4 & 5: Set up and solve the equation for time 't'
# As explained in the plan, the condition for the chain losing contact results in a
# quadratic equation for s = sin(t). The coefficients A, B, and C are derived from
# the geometry and kinematics of the system.
# Equation form: A*s^2 + B*s + C = 0
A = r * sin_alpha
B = -(R * sin_alpha + r * cos_alpha)
C = l_c - h_center - h * cos_alpha - r * sin_alpha

# Step 6: Solve for 's' and then 't'
discriminant = B**2 - 4 * A * C

if discriminant < 0:
    print("No real solution for time exists, which is physically impossible in this scenario.")
else:
    # Solve the quadratic formula for s = sin(t)
    s1 = (-B + math.sqrt(discriminant)) / (2 * A)
    s2 = (-B - math.sqrt(discriminant)) / (2 * A)

    # A valid solution for sin(t) must be in the range [-1, 1].
    valid_s_values = []
    if -1 <= s1 <= 1:
        valid_s_values.append(s1)
    if -1 <= s2 <= 1:
        valid_s_values.append(s2)

    if not valid_s_values:
        print("No valid solution for sin(t) found in the range [-1, 1].")
    else:
        # We need the first time t > 0, which corresponds to the smallest positive sin(t).
        smallest_positive_s = min(s for s in valid_s_values if s >= 0)
        
        # Calculate t from s = sin(t)
        final_time = math.asin(smallest_positive_s)

        # Step 7: Output the results as required
        print("The time 't' is found by solving the following quadratic equation for s = sin(t):")
        print(f"A*s^2 + B*s + C = 0")
        print("\nThe numerical coefficients of the equation are:")
        print(f"A = {A:.4f}")
        print(f"B = {B:.4f}")
        print(f"C = {C:.4f}")
        
        print("\nThus, the final equation to solve is:")
        print(f"{A:.4f} * sin(t)^2 + ({B:.4f}) * sin(t) + {C:.4f} = 0")
        
        print("\nThe first time 't' (in seconds) when the chain loses contact with the ground is:")
        print(final_time)
        
        # The final answer in the required format
        # print(f"\n<<<{final_time:.3f}>>>")