import math

# The problem is to find the limit of the conditional probability
# P[x_0 was not visited | 0 was not visited] for a random walk on a 2D torus.
# Based on the theory of random walks, capacity, and Green's functions,
# the limit can be expressed as 2**(-pi/16).

# The final equation for the limit is: Limit = 2^(-pi/16)
# The numbers in this equation are 2, -1, pi, and 16.

# We calculate the value of this expression.
pi_val = math.pi
exponent = -pi_val / 16
limit_value = 2**exponent

# Output the equation with the numbers and the final computed value.
print("The final equation is: Limit = 2**(-pi / 16)")
print(f"The numbers in the equation are: 2, -1, {pi_val}, 16")
print(f"The value of the limit is: {limit_value}")
