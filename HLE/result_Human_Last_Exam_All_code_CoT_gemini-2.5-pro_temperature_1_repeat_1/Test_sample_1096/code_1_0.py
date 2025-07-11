import math

# We consider a function that is M-smooth and mu-strongly convex.
# The conditioning number is kappa = M/mu.
# The problem asks for the optimal two step sizes gamma_1, gamma_2 for gradient descent
# to minimize the worst-case error ||x_2 - x_*||^2 / ||x_0 - x_*||^2.
# We assume mu=1 and M=kappa, as stated in the problem.

# Let's choose a value for kappa for the calculation.
# A typical value for ill-conditioned problems could be kappa >= 10.
kappa = 10.0

# The optimal step sizes are derived from the roots of the scaled Chebyshev polynomial
# of degree 2. The derivation leads to the following formulas for gamma_1 and gamma_2.
# The order of the step sizes does not matter for the final error, but they are two distinct values.

print(f"For a conditioning number kappa = {kappa}:")
print("\nThe best choice for the pair (gamma_1, gamma_2) is given by the formulas:")
print("gamma_1 = (4*(kappa+1) - 2*sqrt(2)*(kappa-1)) / (kappa^2 + 6*kappa + 1)")
print("gamma_2 = (4*(kappa+1) + 2*sqrt(2)*(kappa-1)) / (kappa^2 + 6*kappa + 1)")
print("\nPlugging in the value of kappa:")

# Denominator of the expression for the step sizes
denominator = kappa**2 + 6*kappa + 1

# Numerators for the two step sizes
# First part of the numerator
num_part1 = 4 * (kappa + 1)
# Second part of the numerator, which is added or subtracted
num_part2 = 2 * math.sqrt(2) * (kappa - 1)

# Calculate the two step sizes
gamma1 = (num_part1 - num_part2) / denominator
gamma2 = (num_part1 + num_part2) / denominator

# Print the equation with numbers for gamma_1
print("\n--- gamma_1 ---")
print(f"gamma_1 = (4*({kappa}+1) - 2*sqrt(2)*({kappa}-1)) / ({kappa}^2 + 6*{kappa} + 1)")
print(f"gamma_1 = ({num_part1} - 2*sqrt(2)*{kappa-1}) / {denominator}")
print(f"gamma_1 = ({num_part1} - {num_part2:.4f}) / {denominator}")
print(f"gamma_1 = {num_part1 - num_part2:.4f} / {denominator}")
print(f"gamma_1 = {gamma1:.6f}")

# Print the equation with numbers for gamma_2
print("\n--- gamma_2 ---")
print(f"gamma_2 = (4*({kappa}+1) + 2*sqrt(2)*({kappa}-1)) / ({kappa}^2 + 6*{kappa} + 1)")
print(f"gamma_2 = ({num_part1} + 2*sqrt(2)*{kappa-1}) / {denominator}")
print(f"gamma_2 = ({num_part1} + {num_part2:.4f}) / {denominator}")
print(f"gamma_2 = {num_part1 + num_part2:.4f} / {denominator}")
print(f"gamma_2 = {gamma2:.6f}")

print("\nNote: The term S mentioned in the prompt, defined as sqrt(M^2 + (M-1)^2), does not appear in the standard derivation for this problem and is likely unrelated.")
