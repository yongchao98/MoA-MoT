import math

# Step 1: Define all given physical parameters.
m = 0.1  # kg (100 g)
M = 10.0 # kg
h = 2.0  # m
theta_deg = 30.0 # degrees
mu = 0.5 # coefficient of friction (1/2)
g = 10.0 # m/s^2

# The initial velocity of the wedge v0 is not needed as it does not affect
# the relative acceleration or the time for relative displacement.

# Step 2: Calculate the relative acceleration 'a_rel' of the block down the incline.
# Convert angle to radians for Python's trigonometric functions
theta_rad = math.radians(theta_deg)
sin_theta = math.sin(theta_rad)
cos_theta = math.cos(theta_rad)

# The formula for relative acceleration is derived from Newton's second law
# for both the block and the wedge.
# a_rel = (g*(sin(θ) - μ*cos(θ))*(M+m)) / (M + m*sin²(θ) + μ*m*sin(θ)*cos(θ))

# Calculate the numerator of the a_rel expression
a_rel_num = g * (sin_theta - mu * cos_theta) * (M + m)

# Calculate the denominator of the a_rel expression
a_rel_den = M + m * sin_theta**2 + m * mu * sin_theta * cos_theta

# Check if the block will slide down (i.e., if the driving force exceeds friction)
if a_rel_num <= 0:
    print("The block does not slide down under these conditions.")
    a_rel = 0.0
else:
    a_rel = a_rel_num / a_rel_den

# Step 3: Calculate the distance 's' the block travels along the incline.
# This is given by s = h / sin(θ).
s = h / sin_theta

# Step 4: Use the kinematic equation to find the time 't'.
# The block starts from rest relative to the wedge, so the equation is:
# s = (1/2) * a_rel * t^2
# Solving for t gives: t = sqrt(2 * s / a_rel)

if a_rel > 0:
    t = math.sqrt(2 * s / a_rel)

    # Step 5: Print the results, including the numbers in the final equation.
    print("To find the time, we use the kinematic equation: t = sqrt(2 * s / a_rel)")
    print("First, we calculate the required values:")
    print(f"Distance to slide (s) = {s:.4f} m")
    print(f"Relative acceleration (a_rel) = {a_rel:.4f} m/s^2")
    
    print("\nPlugging these into the equation:")
    print(f"t = sqrt(2 * {s:.4f} / {a_rel:.4f})")
    
    print("\nThe final result is:")
    print(f"Time = {t:.4f} s")
else:
    # This case will be triggered if a_rel was calculated to be zero or negative.
    print("Time is infinite as the block does not accelerate down the slope.")
