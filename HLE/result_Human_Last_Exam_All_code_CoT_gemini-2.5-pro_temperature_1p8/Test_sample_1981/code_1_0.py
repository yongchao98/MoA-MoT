import decimal

# Set the precision for decimal calculations to ensure accuracy for large numbers.
decimal.getcontext().prec = 100

# From the analysis, we found the value of X0.
# X0^(15/2) = 10^117
# X0 = 10^(117 * 2 / 15) = 10^(78/5)
exponent_X0 = decimal.Decimal('78') / decimal.Decimal('5')
X0 = decimal.Decimal(10) ** exponent_X0

# The expression to be calculated is 10^30 * X0^2 - 10^30 * X0 + 10
thirty = decimal.Decimal(30)
ten = decimal.Decimal(10)

# Calculate each term of the expression
term1 = ten**thirty * X0**2
term2 = ten**thirty * X0
constant_term = ten

# Calculate the final result
result = term1 - term2 + constant_term

# Print the components of the final equation as requested.
# Using string formatting with 'e' for scientific notation.
print(f"{term1:.8e} - {term2:.8e} + {constant_term} = {result:.8e}")

# The problem asks for the final answer in a specific format.
# We will format the final result for the answer block.
final_answer_str = f"{result:.8e}"