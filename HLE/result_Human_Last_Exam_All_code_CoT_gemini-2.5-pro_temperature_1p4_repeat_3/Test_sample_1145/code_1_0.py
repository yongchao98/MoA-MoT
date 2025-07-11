import math

# Define the parameters for the problem
n = 2
x1_ratio = 0.495  # Lower bound of the interval as a fraction of a
x2_ratio = 0.505  # Upper bound of the interval as a fraction of a

# The probability P is the definite integral of |ψ(x)|² = (2/a) * sin²(nπx/a)
# from x1 to x2.
# The analytical result of this integral is [x/a - sin(2nπx/a)/(2nπ)].
# We evaluate this from x/a = x1_ratio to x/a = x2_ratio.

# Let's calculate the components of the final equation
# For n=2, the term in the denominator is 2*n*π = 4π
denominator_term = 2 * n * math.pi

# Value of the antiderivative at the upper bound x2
val_at_x2 = x2_ratio - math.sin(2 * n * math.pi * x2_ratio) / denominator_term

# Value of the antiderivative at the lower bound x1
val_at_x1 = x1_ratio - math.sin(2 * n * math.pi * x1_ratio) / denominator_term

# The probability is the difference
probability = val_at_x2 - val_at_x1

# --- Output the results step-by-step ---
print(f"The probability P is found by solving the definite integral of |ψ(x)|² from {x1_ratio}a to {x2_ratio}a for n={n}.")
print(f"The analytical solution to the integral is P = [x/a - sin({2*n}πx/a)/({2*n}π)] evaluated from {x1_ratio} to {x2_ratio}.")
print("\n--- Plugging in the numbers into the equation ---")

print(f"P = ( {x2_ratio} - sin({2*n}π * {x2_ratio}) / ({denominator_term:.5f}) ) - ( {x1_ratio} - sin({2*n}π * {x1_ratio}) / ({denominator_term:.5f}) )")
print(f"P = ( {val_at_x2:.8f} ) - ( {val_at_x1:.8f} )")
print(f"P = {probability:.8f}")

print(f"\nThe probability that the particle is located in the interval from {x1_ratio}a to {x2_ratio}a for n={n} is approximately {probability:.4e}.")

# Final answer in the required format
# print(f"<<<{probability}>>>")