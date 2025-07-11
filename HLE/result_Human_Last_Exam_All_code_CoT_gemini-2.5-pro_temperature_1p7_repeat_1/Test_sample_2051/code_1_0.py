# Parameters for the specific case in the problem
g = 0
n_plus = 3
n_minus = 1

# The total number of boundaries 'n' is the sum of positively and negatively oriented boundaries.
n = n_plus + n_minus

# The degree of the volume polynomial Z_{g,n} is based on the dimension of the moduli space M_{g,n}.
# The dimension of M_{g,n} is d = 3g - 3 + n.
# The volume is a polynomial in L_i^2 of degree d, which means it is a polynomial
# in L_i of degree 2*d.
# The formula for the degree is: 2 * (3g - 3 + n)

# We will now calculate this value for the given parameters.
print("--- Solution ---")

# (a) Answer regarding continuity
print("(a) Yes. The function Z represents a geometric volume, which varies continuously with its parameters. The polynomial pieces must therefore join continuously.")

# (b) Calculation of the polynomial degree
print("\n(b) The degree of the polynomial Z_{g,n_+,n_-} is calculated as follows:")
print(f"For g = {g}, n+ = {n_plus}, n- = {n_minus}, the total number of boundaries is n = {n_plus} + {n_minus} = {n}.")
print("The degree is given by the formula: 2 * (3*g - 3 + n).")
print("\nSubstituting the values into the formula:")
print(f"Degree = 2 * (3 * {g} - 3 + {n})")

# Calculate the term inside the parentheses
term_in_parentheses = 3 * g - 3 + n
print(f"Degree = 2 * ({term_in_parentheses})")

# Calculate the final degree
degree = 2 * term_in_parentheses
print(f"Degree = {degree}")