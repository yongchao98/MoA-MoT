import decimal

# Set the precision for decimal arithmetic to handle large numbers accurately.
decimal.getcontext().prec = 100

# The derived value for X0 is 10**15.6
X0 = decimal.Decimal(10) ** decimal.Decimal('15.6')

# The expression to calculate is 10**30 * X0**2 - 10**30 * X0 + 10
coeff = decimal.Decimal(10) ** 30
const = decimal.Decimal(10)

# Calculate each term of the equation
term1 = coeff * (X0 ** 2)
term2 = coeff * X0

# Calculate the final result
result = term1 - term2 + const

# As requested, output each number in the final equation
# The numbers are very large, so they are printed in scientific notation by default.
print(f"{term1} - {term2} + {const} = {result}")