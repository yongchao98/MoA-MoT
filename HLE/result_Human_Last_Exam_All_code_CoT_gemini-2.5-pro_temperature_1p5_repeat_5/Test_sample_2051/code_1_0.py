# Parameters for the specific case
g = 0
n_plus = 3
n_minus = 1

# (a) Does the property of piecewise polynomiality imply continuity?
# In the context of volumes of moduli spaces, the function Z is continuous.
# While a general piecewise polynomial function is not always continuous, these
# specific geometric volumes are. The polynomial pieces are defined on cells
# of the moduli space and they join continuously at the boundaries between cells.
answer_a = "Yes"

# (b) Determine the degree of the polynomial Z_{0,3,1}(L_+ | L_-).
# First, calculate the total number of boundaries, n.
n = n_plus + n_minus

# The degree of the polynomial Z_{g,n} is given by the formula: 6g - 6 + 2n.
# We substitute the values g = 0 and n = 4 into this formula.
degree = 6 * g - 6 + 2 * n

# Print the answers.
print("(a) " + answer_a)
print("(b) The degree of the polynomial is found using the formula 6g - 6 + 2n.")
print(f"Given g = {g}, n_+ = {n_plus}, and n_- = {n_minus}, the total number of boundaries is n = {n_plus} + {n_minus} = {n}.")
print("Substituting g and n into the formula:")
print(f"Degree = 6*({g}) - 6 + 2*({n}) = {6*g} - 6 + {2*n} = {degree}")

print(f"<<<(a) {answer_a}; (b) {degree}>>>")