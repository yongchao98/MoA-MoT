# Part (a): Answering the question about continuity.
# While piecewise polynomiality does not universally imply continuity, in the context
# of the volume of moduli spaces Z_{g, n_+, n_-}, the function is known to be continuous.
# The polynomial pieces are guaranteed to meet continuously at the boundaries (walls) of the cells.
print("(a) Yes")
print("-" * 25)

# Part (b): Calculating the degree of the polynomial Z_{0,3,1}.
# The degree of the volume polynomial Z_{g, n_+, n_-} is given by the formula 6g - 6 + 2n,
# where g is the genus and n is the total number of boundaries (n = n_+ + n_-).

# Define the given parameters
g = 0
n_plus = 3
n_minus = 1

# Calculate the total number of boundaries
n = n_plus + n_minus

# Calculate the degree using the formula
degree = 6 * g - 6 + 2 * n

print("(b) To determine the degree of the polynomial Z_{0,3,1}, we use the formula:")
print("Degree = 6g - 6 + 2n, where n = n_plus + n_minus.")
print(f"\nGiven values are:")
print(f"g = {g}")
print(f"n_+ = {n_plus}")
print(f"n_- = {n_minus}")

print(f"\nFirst, calculate the total number of boundaries, n:")
print(f"n = {n_plus} + {n_minus} = {n}")

print(f"\nNow, substitute the values into the degree formula:")
print(f"Degree = 6 * {g} - 6 + 2 * {n}")
print(f"       = {6 * g} - 6 + {2 * n}")
print(f"       = {-6 + 8}")
print(f"       = {degree}")

print(f"\nThe degree of the polynomial is {degree}.")