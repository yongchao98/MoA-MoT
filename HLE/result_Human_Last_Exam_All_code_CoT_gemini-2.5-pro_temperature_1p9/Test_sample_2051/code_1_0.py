# Define parameters for the problem
g = 0
n_plus = 3
n_minus = 1

# (a) Determine if the piecewise polynomial property implies continuity for Z.
# In this specific context of geometric volumes of moduli spaces, the function is continuous.
# The polynomial pieces must join continuously at the cell boundaries.
answer_a = "Yes"

# (b) Determine the degree of the polynomial Z_{0,3,1}.
# The degree of the polynomial Z_{g, n_+, n_-} is given by the formula:
# degree = 4g - 4 + n_+ + n_-.
# We will now calculate this for the given parameters.
degree = 4 * g - 4 + n_plus + n_minus

# Print the final answers as requested.
print(f"(a) {answer_a}")
print(f"(b) The degree is given by the formula 4g - 4 + n_+ + n_-. For g={g}, n_+={n_plus}, n_-={n_minus}, the calculation is:")
# The final part of the request is to "output each number in the final equation"
print(f"Degree = 4*{g} - 4 + {n_plus} + {n_minus} = {degree}")
