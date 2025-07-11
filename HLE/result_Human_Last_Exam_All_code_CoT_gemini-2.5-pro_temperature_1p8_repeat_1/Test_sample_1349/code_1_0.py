import math

# Step 1: Define the components of the final formula for the supremum of X.
# The derivation shows that the supremum of X is given by the formula:
# X_sup = 40 / (36 * pi^2 + 25)

# The numbers in the equation are:
numerator = 40
pi_squared_coefficient = 36
constant_term = 25

# Step 2: Display the formula used for the calculation.
# We are asked to output each number in the final equation.
print("The final expression for the supremum of X is:")
print(f"{numerator} / ({pi_squared_coefficient} * \u03c0\u00b2 + {constant_term})")
print(f"The numbers involved are: numerator = {numerator}, coefficient of \u03c0\u00b2 = {pi_squared_coefficient}, constant = {constant_term}")

# Step 3: Calculate the numerical value.
pi = math.pi
supremum_X = numerator / (pi_squared_coefficient * pi**2 + constant_term)

# Step 4: Print the final numerical answer.
print("\nThe calculated supremum value is:")
print(supremum_X)