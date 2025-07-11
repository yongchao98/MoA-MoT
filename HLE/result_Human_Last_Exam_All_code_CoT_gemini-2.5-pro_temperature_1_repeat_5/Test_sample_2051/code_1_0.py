# Part (a) is a theoretical question. The answer is based on the properties of
# Weil-Petersson volumes of moduli spaces, which are known to be continuous.

# For Part (b), we calculate the degree of the polynomial Z_{0,3,1}.
# The degree of the volume polynomial Z_{g,n} corresponds to the real dimension
# of the moduli space of bordered Riemann surfaces, which is given by the formula 6g - 6 + 2n.

# We will now calculate this for the given parameters.

# Parameters from the problem
g = 0
n_plus = 3
n_minus = 1

# First, calculate the total number of boundaries, n.
n = n_plus + n_minus

# Now, apply the dimension formula to find the degree of the polynomial.
degree = 6 * g - 6 + 2 * n

print("Answer to part (a):")
print("Yes, the property of piecewise polynomiality of Z_{g, n_+, n_-} implies that it is continuous for all L_+ and L_- satisfying the residue condition.")
print("\nCalculation for part (b):")
print(f"The degree of the polynomial Z_{{g,n}} is given by the formula: 6g - 6 + 2n.")
print(f"For Z_{{0,3,1}}, we have g = {g}, n_+ = {n_plus}, and n_- = {n_minus}.")
print(f"The total number of boundaries is n = n_+ + n_- = {n_plus} + {n_minus} = {n}.")
print("\nSubstituting these values into the formula:")
# Output the equation with each number explicitly shown, as requested.
print(f"Degree = (6 * {g}) - 6 + (2 * {n})")
print(f"Degree = {6*g} - 6 + {2*n}")
print(f"Degree = {degree}")
print(f"\nTherefore, the degree of the polynomial Z_{{0,3,1}}(L_+ | L_-) is {degree}.")
