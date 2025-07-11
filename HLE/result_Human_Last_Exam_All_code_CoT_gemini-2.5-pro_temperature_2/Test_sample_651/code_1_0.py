import math

# Based on the derivation, the limit of M(theta) as theta approaches 0 is pi.
# The formula for M(theta) for small theta is M(theta) = pi - theta / 2.
# Let's compute this for a very small theta to demonstrate the limit.
theta = 1e-12
M_theta = math.pi - theta / 2

# The final equation is limit = pi. We will print the numerical value of pi.
final_limit = math.pi
print("The final equation is: M(0) = pi")
print(f"The value is approximately: {final_limit}")
