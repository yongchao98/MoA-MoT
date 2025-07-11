# The problem parameters.
n = 99

# As shown by the mathematical derivation, the sum simplifies to (n+1)^(n-1).
base = n + 1
exponent = n - 1

# The question requires the answer as a power of 10.
# The result is 100^98 = (10^2)^98 = 10^(2 * 98) = 10^196.
final_base = 10
final_exponent = 2 * exponent

# Calculate the final value for confirmation.
final_value = pow(final_base, final_exponent)

# The final equation is 10^196. As requested, we output the numbers
# in this final equation.
print(f"The final calculation leads to the equation: {final_base}^{final_exponent}")
print(f"Base: {final_base}")
print(f"Exponent: {final_exponent}")
print(f"Value: {final_value}")
