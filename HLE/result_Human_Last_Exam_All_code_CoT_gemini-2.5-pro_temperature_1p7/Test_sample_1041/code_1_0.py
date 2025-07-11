import math

# Step 1: Define the parameters of the problem.
R = 1000.0  # Radius of the disk
S = (0.0, 300.0)    # Starting point
a1 = (0.0, 0.0)     # First target point
a2 = (2.0, 0.0)     # Second target point

# The effective radius of a single point in the continuous approximation.
r0 = 1.0

# Step 2: Calculate the effective distance from the start point S to the target set A.
# This is the geometric mean of the distances to the individual points in A.
dist_S_a1 = math.sqrt((S[0] - a1[0])**2 + (S[1] - a1[1])**2)
dist_S_a2 = math.sqrt((S[0] - a2[0])**2 + (S[1] - a2[1])**2)
d_S_A = math.sqrt(dist_S_a1 * dist_S_a2)

# Step 3: Calculate the effective radius of the target set A.
dist_a1_a2 = math.sqrt((a1[0] - a2[0])**2 + (a1[1] - a2[1])**2)
r_eff_A = math.sqrt(r0 * dist_a1_a2)

# Step 4: Apply the formula to compute the probability.
# The formula is P = log(R / d_S_A) / log(R / r_eff_A).
numerator = math.log(R / d_S_A)
denominator = math.log(R / r_eff_A)
probability = numerator / denominator

# Output the components of the final equation and the result.
print(f"The probability is calculated by the equation: log(R/d_S_A) / log(R/r_eff_A)")
print(f"log({R}/{d_S_A}) / log({R}/{r_eff_A})")
print(f"= {numerator} / {denominator}")
# The result to be printed with required precision.
print(f"The final probability is {probability:.3g}")