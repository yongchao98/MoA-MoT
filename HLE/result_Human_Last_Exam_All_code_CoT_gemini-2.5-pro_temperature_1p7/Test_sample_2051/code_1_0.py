# Define the parameters for the specific case
g = 0      # Genus
n_plus = 3 # Number of positively oriented boundaries
n_minus = 1# Number of negatively oriented boundaries

# Calculate the total number of boundaries
n = n_plus + n_minus

# The formula for the degree of the polynomial Z_{g,n} is 6*g - 6 + 2*n.
# We will now calculate this value.
degree = 6 * g - 6 + 2 * n

# Output the explanation and the final result.
# The user wants to see each number in the final equation.
print(f"(a) The property of being a piecewise polynomial volume of a moduli space implies continuity. So the answer is Yes.")
print(f"(b) The degree of the polynomial Z_{{g, n_plus, n_minus}} is calculated using the formula d = 6g - 6 + 2n, where n = n_plus + n_minus.")
print(f"Given g = {g}, n_plus = {n_plus}, and n_minus = {n_minus}, the total number of boundaries is n = {n}.")
print(f"Plugging these values into the formula:")
print(f"Degree = 6 * {g} - 6 + 2 * {n} = {6*g} - 6 + {2*n} = {degree}")
print("\nFinal Answer:")
print(f"(a) Yes; (b) {degree}")