import math

# Step 1: Define constants and determine geometric parameters from the problem statement.
h = 1.0       # m, robot height
r = 0.25      # m, arm length, converted from 25 cm
l_c = 10.0    # m, chain length
d = 20.0      # m, apparent diameter of the path from a side view
l_shadow = 10 * math.sqrt(3) # m, shadow length on the ground

# The radius of the circular path is half its diameter.
R = d / 2.0

# The path is a circle tilted at an angle alpha. With the sun overhead,
# the shadow length l_shadow corresponds to the minor axis of the elliptical shadow.
# l_shadow = 2 * R * cos(alpha)
cos_alpha = l_shadow / (2 * R)
alpha = math.acos(cos_alpha)
sin_alpha = math.sin(alpha)

# Step 2 & 3: Formulate the quadratic equation for x = sin(t).
# The robot's angular speed on the path is omega = v/R = 10/10 = 1 rad/s, so its angle is t.
# The arm's angular speed is dot_beta = 1 rad/s, so its angle is also t.
# The height of the chain's end is z_chain_end(t). Setting it to 0 gives:
# z_chain_end(t) = R*sin(alpha)*(1+sin(t)) + h*cos(alpha) + r*(cos^2(t)*sin(alpha) + sin(t)*cos(alpha)) - l_c = 0
# By substituting cos^2(t) = 1 - sin^2(t) and letting x = sin(t), we get a quadratic equation:
# A*x^2 + B*x + C = 0
# The coefficients are defined as follows:

A = -r * sin_alpha
B = R * sin_alpha + r * cos_alpha
C = R * sin_alpha + h * cos_alpha + r * sin_alpha - l_c

print("The problem is solved by finding the smallest positive time 't' that satisfies z(t)=0.")
print("This can be expressed as a quadratic equation for x = sin(t) of the form A*x^2 + B*x + C = 0.")
print("\nThe numerical coefficients of the equation are:")
print(f"A = {A:.4f}")
print(f"B = {B:.4f}")
print(f"C = {C:.4f}")
print("\nThus, the final equation to solve for sin(t) is:")
print(f"({A:.4f}) * sin^2(t) + ({B:.4f}) * sin(t) + ({C:.4f}) = 0")
print("-" * 30)

# Step 4: Solve the quadratic equation for sin(t) and then find the time t.
discriminant = B**2 - 4 * A * C

final_time = 0
solution_found = False
if discriminant >= 0:
    sqrt_discriminant = math.sqrt(discriminant)
    
    # Calculate the two possible solutions for x = sin(t)
    x1 = (-B + sqrt_discriminant) / (2 * A)
    x2 = (-B - sqrt_discriminant) / (2 * A)
    
    # Find the valid physical solution for sin(t), which must be in the range [-1, 1].
    valid_solutions = []
    if -1.0 <= x1 <= 1.0:
        valid_solutions.append(x1)
    if -1.0 <= x2 <= 1.0:
        valid_solutions.append(x2)

    if valid_solutions:
        # We need the smallest positive time t, so we need the smallest positive sin(t) value.
        # Since z(0) < 0 and z'(0) > 0, the first root for t > 0 corresponds to the smallest positive root for sin(t).
        sin_t_sol = min(s for s in valid_solutions if s >= 0)
        
        # Calculate time t = arcsin(sin(t)).
        # For a positive sin(t), math.asin gives the smallest positive time t.
        final_time = math.asin(sin_t_sol)
        print(f"The first time the chain loses contact with the ground is at t = {final_time:.4f} seconds.")
        solution_found = True

if not solution_found:
    print("No physical solution found for when the chain leaves the ground.")

print(f"\n<<< {final_time:.4f} >>>")