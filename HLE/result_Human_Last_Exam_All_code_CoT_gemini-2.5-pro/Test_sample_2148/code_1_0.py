import math

# Step 1: Define the coupling constants for neutrinos
c_V = 0.5
c_A = 0.5

# Step 2: Calculate the common term (c_V^2 + c_A^2)
couplings_sq_sum = c_V**2 + c_A**2

# Step 3: Calculate X1 and X2 based on the derived formulas
# X1 = 32 * sqrt(2) * (c_V^2 + c_A^2)
# X2 = (c_V^2 + c_A^2) / (12 * sqrt(2) * pi)
X1 = 32 * math.sqrt(2) * couplings_sq_sum
X2 = couplings_sq_sum / (12 * math.sqrt(2) * math.pi)

# Step 4: Calculate the product X1 * X2
product = X1 * X2

# Step 5: Calculate the inverse of the product.
# The symbolic result is 3*pi/2
final_result = 1 / product
numerator = 3
denominator = 2

# Step 6: Print the final equation and its numerical value
print("The final equation is: (X1 * X2)^-1 = ({} * pi) / {}".format(numerator, denominator))
print("The numerical value is: {:.15f}".format(final_result))

<<<4.712388980384690>>>