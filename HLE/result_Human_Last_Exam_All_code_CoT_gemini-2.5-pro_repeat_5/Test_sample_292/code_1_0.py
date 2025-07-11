import math

# Define the value of n from the problem statement.
n = 99

# The analytical derivation shows the sum is equal to (n+1)^(n-1).
base = n + 1
exponent = n - 1

# The problem requires the answer as a power of 10.
# We can express the result as follows:
# (n+1)^(n-1) = 100^98 = (10^2)^98 = 10^(2*98) = 10^196
final_base = 10
final_exponent = 2 * exponent

# Print the final equation with each number, as requested.
# This shows how the expression simplifies to a power of 10.
print(f"The sum is equal to ({n}+1)^({n}-1), which simplifies to a power of 10 as shown in the equation below:")
print(f"{base}^{exponent} = ({final_base}^2)^{exponent} = {final_base}^(2 * {exponent}) = {final_base}^{final_exponent}")
