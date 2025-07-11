import math

# Step 1: Define the parameters of the problem based on the plan.

# R is the radius of the disk the random walk is confined to.
R = 1000.0

# The starting point of the random walk.
start_point = (0, 300)

# The target set is {(0,0), (2,0)}. We find its geometric center.
target_center = ((0 + 2) / 2.0, (0 + 0) / 2.0) # which is (1.0, 0.0)

# d is the distance from the starting point to the center of the target set.
d = math.sqrt((start_point[0] - target_center[0])**2 + (start_point[1] - target_center[1])**2)

# The target set is two points at distance L=2 apart.
# We approximate its effective radius (r_eff) using the logarithmic capacity
# of a line segment of length L, which is L/4.
L = 2.0
r_eff = L / 4.0

# Step 2: Apply the hitting probability formula.
# P = ln(R / d) / ln(R / r_eff)
# We use the natural logarithm (math.log). Any logarithm base could be used
# as it cancels out, but natural log is standard.
numerator = math.log(R / d)
denominator = math.log(R / r_eff)
probability = numerator / denominator

# Step 3: Print the details of the calculation and the result.
# The problem requests to output each number in the final equation.
print("The probability P is estimated by the formula: P = log(R / d) / log(R / r_eff)")
print("\nParameters:")
print(f"R (Disk Radius) = {R}")
print(f"d (Distance to Target) = sqrt(({start_point[0] - target_center[0]})**2 + ({start_point[1] - target_center[1]})**2) = {d}")
print(f"r_eff (Effective Target Radius) = {L} / 4 = {r_eff}")

print("\nFinal Equation:")
# Display the equation with the calculated numbers
print(f"P = log({R} / {d:.4f}) / log({R} / {r_eff})")
print(f"P = log({R/d:.4f}) / log({R/r_eff:.4f})")
print(f"P = {numerator:.4f} / {denominator:.4f}")

# Print the final result formatted to three significant digits.
print(f"\nThe probability is approximately: {probability:.3f}")
