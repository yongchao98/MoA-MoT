import math

# Step 1: Define the given constants
T = math.log(10)

# Step 2: Calculate the terms involving T
# e_minus_2T = 10^(-2) = 0.01
# e_minus_3T = 10^(-3) = 0.001
term1 = 1 - math.exp(-2 * T)  # This is 1 - 0.01 = 0.99
term2 = 1 - math.exp(-3 * T)  # This is 1 - 0.001 = 0.999

# Step 3: Calculate alpha cubed based on the derived formula
# alpha^3 = (3 * (1 - exp(-2T))^4) / (16 * (1 - exp(-3T)))
alpha_cubed = (3 * (term1**4)) / (16 * term2)

# Step 4: Calculate alpha by taking the cube root
alpha = alpha_cubed**(1/3)

# Print the final numerical result for alpha
print(f"The calculated value of alpha is: {alpha}")

# As requested by the prompt, output the numbers in the final equation used to determine alpha.
# The equation is alpha^3 = (3 * (1 - e^(-2T))^4) / (16 * (1 - e^(-3T)))
# Here, we print the components of this equation.
print("\nBreakdown of the calculation for alpha:")
print(f"T = {T}")
print(f"The term (1 - e^(-2T)) is: {term1}")
print(f"The term (1 - e^(-3T)) is: {term2}")
print(f"The numerator in the expression for alpha^3 is: 3 * ({term1})^4 = {3 * (term1**4)}")
print(f"The denominator in the expression for alpha^3 is: 16 * {term2} = {16 * term2}")
print(f"alpha^3 = {alpha_cubed}")
print(f"alpha = (alpha^3)^(1/3) = {alpha}")
