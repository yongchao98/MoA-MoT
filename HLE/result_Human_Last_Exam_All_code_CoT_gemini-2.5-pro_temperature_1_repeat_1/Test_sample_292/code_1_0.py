# The problem asks for a sum which, through mathematical derivation,
# simplifies to a clean expression.
# Let n be the vocabulary size, n = 99.
# The sum is equal to (n + 1)^(n - 1).

n = 99

# The final equation for the sum S is S = (n+1)^(n-1)
base = n + 1
exponent = n - 1

# Calculate the result
result = pow(base, exponent)

# The problem asks to write the answer as a power of 10.
# S = 100^98 = (10^2)^98 = 10^(2*98) = 10^196
power_of_10_exponent = 2 * exponent

print(f"The final equation for the sum is: S = {base}^{exponent}")
print(f"The value of the sum is: {result}")
print(f"As a power of 10, the result is: 10^{power_of_10_exponent}")