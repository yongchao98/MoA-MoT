import math

# The diameter of the quarter-sphere package is given.
D_package = 250.0  # in cm

# As derived in the plan, the formula for the maximum diameter (d) of the inscribed sphere is:
# d = D * (sqrt(2) - 1)
# where D is the diameter of the package.

# Calculate the values needed for the equation.
sqrt_2 = math.sqrt(2)
term_in_parentheses = sqrt_2 - 1
max_diameter = D_package * term_in_parentheses

# Print the explanation and the steps of the calculation as requested.
print("The formula for the maximum diameter 'd' of a sphere that fits inside a quarter-sphere package of diameter 'D' is:")
print("d = D * (sqrt(2) - 1)")
print("\n--- Calculation Steps ---")
print(f"1. The given diameter of the package is D = {D_package} cm.")
print(f"2. The value of sqrt(2) is approximately {sqrt_2}.")
print(f"3. Plugging the values into the formula:")
print(f"   d = {D_package} * ({sqrt_2} - 1)")
print(f"4. First, evaluate the expression in the parentheses:")
print(f"   d = {D_package} * {term_in_parentheses}")
print(f"5. Finally, perform the multiplication to find the diameter:")
print(f"   d = {max_diameter}")
print("\n-------------------------")
print(f"The maximum diameter of the sphere is {max_diameter} cm.")
