import math

# Step 1: Define the equation for the generating amplitudes.
# Based on the analysis of the van der Pol equation, the equation connecting
# the amplitudes c1 and c2 of the generating solution is c1^2 + c2^2 = 4.
# This equation describes the amplitude of the stable limit cycle.
equation_str = "c1^2 + c2^2 = 4"
print(f"The equation for the generating amplitudes is: {equation_str}")

# Step 2: Apply the condition c1 = c2.
# Substituting c2 with c1 in the equation gives:
# c1^2 + c1^2 = 4
# 2 * c1^2 = 4
# This is the final equation we need to solve for c1.
a = 2
b = 4
final_equation_str = f"{a} * c1^2 = {b}"
print(f"Given the condition c1 = c2, the final equation to solve is: {final_equation_str}")
print(f"The numbers in this final equation are {a} and {b}.")

# Step 3: Solve the equation for c1.
# c1^2 = 4 / 2 = 2
c1_squared = b / a

# We need the first positive root, so we take the positive square root.
c1 = math.sqrt(c1_squared)

# Step 4: Print the final result.
print(f"The positive solution for c1 is sqrt({c1_squared}), which is {c1}.")