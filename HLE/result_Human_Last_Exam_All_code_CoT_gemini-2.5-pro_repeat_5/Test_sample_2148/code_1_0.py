import math

# Step 1: Define the coupling constants for neutrinos
c_V = 0.5
c_A = 0.5

# Calculate the sum of squares of the couplings
c_V2_plus_c_A2 = c_V**2 + c_A**2

# Step 2: Calculate X2 based on the formulas provided in the problem
# X2 = (c_V^2 + c_A^2) / (12 * sqrt(2) * pi)
X2 = c_V2_plus_c_A2 / (12 * math.sqrt(2) * math.pi)

# Step 3: Calculate X1 using the derived physical relationship X1 = 16 * pi * X2
X1 = 16 * math.pi * X2

# Step 4: Calculate the final result (X1 * X2)^-1
# The final equation is R = (X1 * X2)^-1
# The numbers in this equation are X1, X2, and the result R.

product = X1 * X2
result = 1 / product

print(f"The value of X1 is: {X1}")
print(f"The value of X2 is: {X2}")
print(f"The product (X1 * X2) is: {product}")
print(f"The final result for (X1 * X2)^-1 is: {result}")
print("\nFor comparison, the symbolic result is 72 * pi.")
print(f"The numerical value of 72 * pi is: {72 * math.pi}")
