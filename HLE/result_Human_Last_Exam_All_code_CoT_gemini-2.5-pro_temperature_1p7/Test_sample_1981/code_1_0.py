import decimal

# Set a high precision for the decimal calculations to handle the large numbers
decimal.getcontext().prec = 100

# From our derivation, X_0 = 10^(78/5)
# We use the Decimal type for accurate floating-point arithmetic.
ten = decimal.Decimal(10)
thirty = decimal.Decimal(30)
seventy_eight = decimal.Decimal('78')
five = decimal.Decimal('5')

# Calculate X_0
X0_exponent = seventy_eight / five
X0 = ten**X0_exponent

# Calculate X_0^2
X0_sq = X0**2

# Define the coefficients in the final expression
c1 = ten**thirty
c2 = -ten**thirty
c3 = ten

# Calculate the final result
result = c1 * X0_sq + c2 * X0 + c3

# As requested, we print out each number in the final equation.
# Since the numbers are extremely large, we format them in scientific notation.
print(f"The expression to evaluate is: ({c1:e}) * ({X0_sq:e}) + ({c2:e}) * ({X0:e}) + ({c3})")
print(f"This evaluates to: {result:e}")