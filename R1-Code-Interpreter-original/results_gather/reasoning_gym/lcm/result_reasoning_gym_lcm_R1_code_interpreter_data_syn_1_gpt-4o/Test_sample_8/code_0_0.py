import math

# Define the numbers
a = 17
b = 98

# Calculate the GCD
gcd = math.gcd(a, b)

# Calculate the LCM using the formula
lcm = abs(a * b) // gcd

# Print the result
print(lcm)