import math

# Define the numbers
a = 33
b = 73

# Calculate the GCD
gcd = math.gcd(a, b)

# Calculate the LCM using the relationship between GCD and LCM
lcm = abs(a * b) // gcd

# Output the result
print(lcm)