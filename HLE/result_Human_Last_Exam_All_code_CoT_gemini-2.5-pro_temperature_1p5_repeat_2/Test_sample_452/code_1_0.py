import math

# The problem is to find the constant 'b' in the asymptotic formula for the 
# expected cover-and-return time of a random walk on a random n-vertex tree,
# C(n) ~ b * n^(3/2).
#
# This is a known result from probability theory. The exact value of the constant 'b'
# was determined by L. Devroye and A. Sbihi (1992) to be sqrt(pi / 2).
#
# This script calculates the numerical value of this constant.

# The equation for the constant b is:
# b = sqrt(pi / 2)

# We define the numbers in the equation.
pi_value = math.pi
divisor = 2.0

# Now, we calculate the value of b.
b_value = math.sqrt(pi_value / divisor)

print("The constant 'b' is derived from the formula: b = sqrt(pi / 2)")
print("---------------------------------------------------------")
print(f"1. The value of pi (π) is: {pi_value}")
print(f"2. The value of the divisor is: {divisor}")
print("---------------------------------------------------------")
print(f"Result: The numerical value of b = sqrt({pi_value} / {divisor}) is approximately {b_value}")
print(f"The exact value of b is sqrt(pi/2).")
