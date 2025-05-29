import math

# Define the numbers
a = 48
b = 99

# Calculate the GCD
gcd = math.gcd(a, b)

# Calculate the LCM using the relationship between LCM and GCD
lcm = abs(a * b) // gcd

# Print the result
print(lcm)