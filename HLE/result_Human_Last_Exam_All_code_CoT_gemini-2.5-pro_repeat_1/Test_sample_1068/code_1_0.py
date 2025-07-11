import math

# The value of p provided in the problem.
p = 14

# Based on the step-by-step derivation, the complex integral simplifies to the formula:
# l(p) = 2 * p * ln(2)

# We can now calculate the value for p = 14.
# Let's define the numbers in the final equation.
factor1 = 2
factor2 = p
ln_argument = 2

# Calculate the coefficient and the logarithm
coefficient = factor1 * factor2
log_value = math.log(ln_argument)

# Calculate the final result
result = coefficient * log_value

# Output the final equation with each number explicitly shown, as requested.
print("The final simplified equation is l(p) = 2 * p * ln(2).")
print(f"For p = {p}, the equation with its numbers is:")
print(f"l({p}) = {factor1} * {p} * ln({ln_argument})")
print(f"l({p}) = {coefficient} * ln({ln_argument})")
print("\nThe final numerical value is:")
print(result)
