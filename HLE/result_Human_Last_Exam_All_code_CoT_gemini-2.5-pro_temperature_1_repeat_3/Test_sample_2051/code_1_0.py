# This script answers the two-part question about the volume of moduli spaces of ribbon graphs.

# Part (a): Analysis of continuity.
# While in general a piecewise polynomial function is not guaranteed to be continuous,
# the function Z in this context represents a geometric volume. Volumes are inherently
# continuous with respect to their defining parameters. Therefore, the function Z is continuous.
# The piecewise polynomial nature means that the polynomial expressions for different
# regions (cells) must join continuously at their boundaries.
answer_a = "Yes"

# Part (b): Calculation of the polynomial degree.
# The degree of the volume polynomial Z_{g,n} as a function of the boundary lengths L_i
# is given by the formula: 2 * (3g - 3 + n), where g is the genus and n is the total
# number of boundaries.

# We are given the following parameters for Z_{g, n_+, n_-}:
g = 0
n_plus = 3
n_minus = 1

# First, we calculate the total number of boundaries, n.
n = n_plus + n_minus

# Next, we calculate the term inside the parentheses, which corresponds to the
# complex dimension of the moduli space of curves M_{g,n}.
dim_term = 3 * g - 3 + n

# Finally, we calculate the degree of the polynomial.
degree = 2 * dim_term

# Now, we print the results in a clear format.
print(f"(a) {answer_a}")
print(f"(b) The degree of the polynomial Z_{{0,3,1}} is {degree}.")
print("\n--- Calculation for (b) ---")
print(f"The degree is determined by the formula: Degree = 2 * (3*g - 3 + n)")
print(f"Given parameters: g = {g}, n_+ = {n_plus}, n_- = {n_minus}")
print(f"Total number of boundaries n = {n_plus} + {n_minus} = {n}")
print("Substituting the values into the formula:")
print(f"Degree = 2 * (3*{g} - 3 + {n})")
print(f"Degree = 2 * ({dim_term})")
print(f"Degree = {degree}")