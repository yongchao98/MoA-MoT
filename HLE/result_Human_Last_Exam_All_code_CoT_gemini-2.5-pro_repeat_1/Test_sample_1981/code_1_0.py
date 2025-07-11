import math

# Step 1: Calculate X0
# From the derivation, X0^(15/2) = 10^117.
# So, X0 = 10^(117 * 2 / 15) = 10^(15.6).
x0_exp = 15.6
X0 = 10**x0_exp

# Step 2: Calculate the terms of the final expression
# The expression is 10^30 * X0^2 - 10^30 * X0 + 10.
coeff = 10**30
X0_squared = X0**2
const = 10

term1 = coeff * X0_squared
term2 = coeff * X0

# Step 3: Calculate the final result
result = term1 - term2 + const

# As per the instruction "output each number in the final equation!",
# we print the components of the calculation.
# The equation is: 10^30 * (10^15.6)^2 - 10^30 * 10^15.6 + 10 = result
# which is: 10^61.2 - 10^45.6 + 10 = result

print("The equation to solve is: C * X0^2 - C * X0 + D")
print(f"C = 10**30")
print(f"X0 = 10**15.6 = {X0:.4e}")
print(f"D = 10")
print("\nSubstituting the values:")
# Using f-string formatting with exponents for clarity
print(f"10**30 * (10**{x0_exp})**2 - 10**30 * 10**{x0_exp} + 10")
print(f"= 10**{30 + x0_exp*2} - 10**{30 + x0_exp} + 10")
print(f"= 10**{61.2} - 10**{45.6} + 10")
print(f"\nFinal calculated result: {result:.4e}")

# The final result is a large floating-point number.
# Due to the large difference in magnitude, 10**61.2 - 10**45.6 is numerically
# indistinguishable from 10**61.2 in standard floating-point arithmetic.
# The result is effectively 10**61.2 + 10.
final_answer = 10**61.2 + 10
# We will use the result from python's calculation which handles precision.
final_answer_val = result
print(f"Final Answer is {final_answer_val}")
