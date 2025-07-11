from decimal import Decimal, getcontext

# Set precision for Decimal calculations
getcontext().prec = 100

# Step 5: Calculate X_0
# X_0^(15/2) = 10^117
# X_0 = (10^117)^(2/15) = 10^(15.6)
X0 = Decimal(10) ** Decimal('15.6')

# Step 6: Calculate the final expression 10^30 * X_0^2 - 10^30 * X_0 + 10
X0_sq = X0 ** 2
term1 = Decimal(10)**30 * X0_sq
term2 = Decimal(10)**30 * X0
term3 = Decimal(10)

# The "final equation" consists of these numbers
print(f"The equation is: {term1} - {term2} + {term3}")

# Calculate the final result
result = term1 - term2 + term3

print(f"\nThe value of X_0 is: {X0}")
print(f"The value of X_0^2 is: {X0_sq}")
print(f"The final result is: {result}")
print(f"\nFinal Answer in <<<...>>> format:")
print(f"<<<{result}>>>")