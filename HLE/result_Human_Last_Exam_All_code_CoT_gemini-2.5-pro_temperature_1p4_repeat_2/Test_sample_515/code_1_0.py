import math

# The probability is given by the formula 1 / (2 * ln(3600) + gamma + ln(8))
# We need to calculate the value of the denominator.

# The point y is (3600, 0)
n = 3600

# Natural logarithm of n
ln_n = math.log(n)

# Euler-Mascheroni constant
gamma = 0.5772156649

# Natural logarithm of 8
ln_8 = math.log(8)

# Calculate the denominator of the probability formula
denominator = 2 * ln_n + gamma + ln_8

# The probability is the reciprocal of the denominator
probability = 1 / denominator

# Print the components of the calculation
print(f"The starting point is x_0 = (0,1). The target set A is the set of neighbors of y = (3600,0).")
print(f"The probability is approximated by a(x_0)/a(y).")
print(f"a(x_0) = a(0,1) = a(1,0) = 1/pi")
# Using the asymptotic formula for a(y) from Spitzer
print(f"a(y) = a(3600,0) is approximated by (1/pi) * (2*ln(3600) + gamma + ln(8)).")
print(f"The ratio simplifies to 1 / (2*ln(3600) + gamma + ln(8)).")
print("\n--- Calculation ---")
print(f"n = {n}")
print(f"2 * ln(n) = 2 * ln({n}) = {2 * ln_n:.4f}")
print(f"gamma = {gamma:.4f}")
print(f"ln(8) = {ln_8:.4f}")
print(f"Denominator = {2 * ln_n:.4f} + {gamma:.4f} + {ln_8:.4f} = {denominator:.4f}")
print(f"Probability = 1 / {denominator:.4f} = {probability:.4f}")
print(f"\nThe approximate probability is {probability:.2g}.")

final_equation = f"1 / (2 * ln({n}) + {gamma} + {ln(8)})"
final_result = 1/(2*math.log(n)+gamma+math.log(8))
print(f"\nFinal Equation: P = {final_equation} = {final_result}")
