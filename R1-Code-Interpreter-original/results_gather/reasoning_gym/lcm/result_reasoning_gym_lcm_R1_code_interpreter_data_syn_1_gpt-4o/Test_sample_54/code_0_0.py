import math

# Define the numbers
a = 99
b = 9

# Calculate the GCD
gcd = math.gcd(a, b)

# Calculate the LCM using the relationship between GCD and LCM
lcm = abs(a * b) // gcd

# Print the result
print(lcm)