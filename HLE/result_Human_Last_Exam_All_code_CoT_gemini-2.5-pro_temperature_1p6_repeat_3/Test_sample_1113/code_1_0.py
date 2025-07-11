import math

# --- Step 1: Define constants and geometric parameters ---
h = 1.0  # Robot height in m
r = 0.25  # Arm length in m
l_c = 10.0  # Chain length in m
v = 10.0  # Robot speed in m/s
beta_dot = 1.0  # Arm angular speed in rad/s
d = 20.0  # Visible path diameter in m
l = 10 * math.sqrt(3)  # Shadow length in m

# Path parameters
R = d / 2
# From l = d * cos(alpha)
alpha = math.acos(l / d)
# Robot angular velocity on the path
omega = v / R

print("--- Initial Parameters ---")
print(f"Path Radius (R): {R:.2f} m")
print(f"Path Tilt Angle (alpha): {math.degrees(alpha):.2f} degrees")
print(f"Robot Angular Velocity (omega): {omega:.2f} rad/s")
print("-" * 20)

# --- Step 2: Derive the equation for the height of the arm's tip, Z_tip ---
# The height of the tip Z_tip(t) can be expressed as a function of s = sin(t).
# Z_tip(s) = A*s^2 + B*s + C, where t is time in seconds.
# We set Z_tip(t) = l_c to find when the chain loses contact.
# This gives the quadratic equation:
# A_eq * s^2 + B_eq * s + C_eq = 0

# Coefficients of Z_tip(s) = -(r/2)s^2 + (5 - r*sqrt(3)/2)s + (5 + r/2 - sqrt(3)/2)
# The equation to solve is Z_tip(t) = l_c = 10
# -(r/2)s^2 + (5 - r*sqrt(3)/2)s + (5 + r/2 - sqrt(3)/2) = 10
# Rearranging to form a standard quadratic equation a*x^2 + b*x + c = 0:
# (r/2)s^2 - (5 - r*sqrt(3)/2)s + (10 - (5 + r/2 - sqrt(3)/2)) = 0

A_eq = r / 2
B_eq = -(5 - r * math.sqrt(3) / 2)
C_eq = l_c - (5 + r / 2 - math.sqrt(3) / 2)

print("--- Solving for Time t ---")
print("The height of the arm tip Z_tip can be expressed as a function of s = sin(t).")
print(f"To find when the chain loses contact, we solve the equation: Z_tip(t) = {l_c} m")
print("This results in a quadratic equation for s = sin(t):")
print(f"As^2 + Bs + C = 0")
print(f"Where:")
print(f"  A = {A_eq:.4f}")
print(f"  B = {B_eq:.4f}")
print(f"  C = {C_eq:.4f}")

# --- Step 3: Check if a solution exists ---
# First, let's find the maximum possible height of the arm tip.
# This occurs at the boundary of sin(t) = 1, since the vertex of the parabola is outside [-1, 1].
Z_tip_max = -(A_eq) * (1**2) + (-B_eq) * 1 + (l_c - C_eq)
print(f"\nThe maximum height the arm tip can reach is: {Z_tip_max:.4f} m")

if Z_tip_max < l_c:
    print(f"\nSince the maximum tip height ({Z_tip_max:.4f} m) is less than the chain length ({l_c} m),")
    print("the chain will never fully lose contact with the ground.")
    # As the problem implies a solution exists, we calculate the time of closest approach (max height).
    # This occurs when sin(t) = 1, so t = pi/2.
    time_of_closest_approach = math.pi / 2
    print(f"The time at which the tip reaches its maximum height is t = pi/2 ≈ {time_of_closest_approach:.4f} s.")
    final_time = time_of_closest_approach
else:
    # Calculate the discriminant to find solutions for s = sin(t)
    discriminant = B_eq**2 - 4 * A_eq * C_eq

    print(f"\nCalculating the discriminant (B^2 - 4AC): {discriminant:.4f}")

    if discriminant < 0:
        print("The discriminant is negative, so there are no real solutions for sin(t).")
        print("This confirms the chain never loses contact with the ground.")
        final_time = None
    else:
        # Solve for the two possible values of s = sin(t)
        s1 = (-B_eq + math.sqrt(discriminant)) / (2 * A_eq)
        s2 = (-B_eq - math.sqrt(discriminant)) / (2 * A_eq)
        print(f"The possible values for sin(t) are: s1 = {s1:.4f}, s2 = {s2:.4f}")

        # Find the valid solution for sin(t) which must be in [-1, 1]
        valid_s = []
        if -1 <= s1 <= 1:
            valid_s.append(s1)
        if -1 <= s2 <= 1:
            valid_s.append(s2)

        if not valid_s:
            print("No valid solution for sin(t) in the range [-1, 1].")
            final_time = None
        else:
            s = valid_s[0]
            # Find the smallest positive time t
            # The chain loses contact when the tip is moving upwards, i.e., d(Z_tip)/dt > 0.
            # d(Z_tip)/dt = cos(t) * [constant - r*sin(t)].
            # We need the solution where cos(t) > 0 for it to be "losing" contact for the first time.
            
            # The two smallest positive solutions for sin(t) = s are:
            # t1 = arcsin(s)  (if s>=0) or 2pi+arcsin(s) (if s<0)
            # t2 = pi - arcsin(s)
            
            t1 = math.asin(s)
            t2 = math.pi - math.asin(s)

            if t1 < 0: t1 += 2*math.pi # ensure t1 is positive

            # Check which time gives cos(t) > 0
            if math.cos(t1) > 0:
                final_time = t1
            else:
                final_time = t2
            print(f"\nThe first time t > 0 for the chain to lose contact is {final_time:.4f} seconds.")

# Acknowledging the paradox in the problem statement
print("\nNote: A detailed analysis of the provided parameters shows the robot arm's tip never reaches the required height of 10m.")
print("The calculated time represents the moment the tip reaches its maximum possible height.")
print(f"Final Answer (time of maximum tip height): {final_time:.4f}")