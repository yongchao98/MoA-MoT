import math

# Define the parameters of the problem
R = 1000.0
S = (0, 300)
A1 = (0, 0)
A2 = (2, 0)
r0 = 1.0  # Effective radius of a single lattice point

# Calculate the required distances
dist_S_A1 = math.sqrt((S[0] - A1[0])**2 + (S[1] - A1[1])**2)
dist_S_A2 = math.sqrt((S[0] - A2[0])**2 + (S[1] - A2[1])**2)
dist_A1_A2 = math.sqrt((A1[0] - A2[0])**2 + (A1[1] - A2[1])**2)

# Calculate the terms for the numerator of the formula
num_term1 = math.log(R / dist_S_A1)
num_term2 = math.log(R / dist_S_A2)
numerator = num_term1 + num_term2

# Calculate the terms for the denominator of the formula
den_term1 = math.log(R / r0)
den_term2 = math.log(R / dist_A1_A2)
denominator = den_term1 + den_term2

# Calculate the final probability
probability = numerator / denominator

# Print the explanation and the final equation with values
print("The probability is estimated using the continuous approximation for a 2D random walk (potential theory).")
print("The formula for the probability h(S) is:")
print("h(S) = (log(R/|S-A1|) + log(R/|S-A2|)) / (log(R/r0) + log(R/|A1-A2|))\n")

print("Plugging in the values:")
print(f"R = {R}")
print(f"S = {S}")
print(f"A1 = {A1}, A2 = {A2}")
print(f"|S-A1| = {dist_S_A1:.4f}")
print(f"|S-A2| = {dist_S_A2:.4f}")
print(f"|A1-A2| = {dist_A1_A2:.4f}")
print(f"r0 = {r0} (effective radius of a point)\n")

print("The equation with the computed values is:")
print(f"Probability = (log({R}/{dist_S_A1:.4f}) + log({R}/{dist_S_A2:.4f})) / (log({R}/{r0}) + log({R}/{dist_A1_A2:.4f}))")
print(f"Probability = ({num_term1:.5f} + {num_term2:.5f}) / ({den_term1:.5f} + {den_term2:.5f})")
print(f"Probability = {numerator:.5f} / {denominator:.5f}")
print(f"Probability = {probability:.5f}\n")

print(f"The final probability rounded to three significant digits is: {probability:.3f}")
<<<0.183>>>