import math

# Step 1: Define the size of the vocabulary from the problem.
n = 99

# Step 2: Use the simplified formula for the sum: (n+1)^(n-1).
# For n = 99, this is (99+1)^(99-1) = 100^98.
base_100 = n + 1
power_98 = n - 1

# Step 3: Express the result as a power of 10.
# 100^98 = (10^2)^98 = 10^(2 * 98)
base_10 = 10
final_exponent = 2 * power_98

# Step 4: Output the numbers in the final equation as requested.
# This shows the result both as a power of 100 and as a power of 10.
print(f"{base_100}^{power_98} = {base_10}^{final_exponent}")