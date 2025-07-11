import math

# Step 1: Define the parameters based on the problem statement.
# R is the radius of the disk.
R = 1000.0
# r is the distance from the starting point (0, 300) to the center of the target set {(0,0), (2,0)}, which is (1,0).
# The distance is sqrt((0-1)^2 + (300-0)^2) = sqrt(1 + 90000), which is very close to 300.
r = 300.0
# d is the distance between the two points in the target set.
d = 2.0
# epsilon is the effective radius (logarithmic capacity) of the two-point target set, which is sqrt(d).
epsilon = math.sqrt(d)

# Step 2: Calculate the probability using the approximation formula.
# P = ln(R/r) / ln(R/epsilon)
# We use the natural logarithm (ln), as the base of the logarithm does not affect the ratio.
numerator = math.log(R / r)
denominator = math.log(R / epsilon)
probability = numerator / denominator

# Step 3: Print the calculation steps and the final result.
print("The probability P is approximated by the formula: P = ln(R/r) / ln(R/epsilon)")
print(f"Here, R = {R}, r = {r}, and epsilon = sqrt({d}) = {epsilon:.4f}")
print(f"P = ln({R}/{r}) / ln({R}/{epsilon:.4f})")
print(f"P = ln({R/r:.4f}) / ln({R/epsilon:.4f})")
print(f"P = {numerator:.4f} / {denominator:.4f}")
print(f"The calculated probability is: {probability:.5f}")

# Final answer with three significant digits.
final_answer = round(probability, 3)
print(f"The probability rounded to three significant digits is: {final_answer}")
