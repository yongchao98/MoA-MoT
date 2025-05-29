import math

# Define the numbers
a = 26
b = 94

# Calculate the GCD
gcd = math.gcd(a, b)

# Calculate the LCM using the relationship between GCD and LCM
lcm = abs(a * b) // gcd

# Print the LCM
print(lcm)