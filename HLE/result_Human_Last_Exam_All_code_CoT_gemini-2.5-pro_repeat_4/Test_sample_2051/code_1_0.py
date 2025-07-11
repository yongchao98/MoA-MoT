import sympy

# --- Final Answer to Part (a) ---
# The question: Does the property of piecewise polynomiality of Z imply that it is continuous?
# The answer is Yes. In the context of geometric volumes, the function is continuous.
# The different polynomial pieces match at the boundaries of the cells.
print("(a) Yes")
print("-" * 20)

# --- Calculation for Part (b) ---
# The question: For g=0, n_+=3, n_-=1, determine the degree of the polynomial Z.
print("(b) Calculation of the degree:")

# Parameters for the problem
g = 0
n_plus = 3
n_minus = 1
n = n_plus + n_minus

# The degree of the volume polynomial Z_{g,n} in the boundary lengths L_i
# is given by the formula 2 * (3g - 3 + n).
degree_formula = 2 * (3 * g - 3 + n)
print(f"The degree of the base polynomial Z_{{{g},{n}}} is given by the formula 2*(3g-3+n).")
print(f"For g={g} and n={n}, the degree is 2 * (3*{g} - 3 + {n}) = {degree_formula}.")

# Define the symbolic variables for the boundary lengths
L1, L2, L3, L4 = sympy.symbols('L1 L2 L3 L4')

# The polynomial for Z_{0,4} is known to be Z = 1/2 * (L1^2 + L2^2 + L3^2 + L4^2)
Z_0_4 = (L1**2 + L2**2 + L3**2 + L4**2) / 2
print(f"\nThe polynomial Z_{{0,4}} is: {Z_0_4}")

# The residue condition for (n_+, n_-) = (3,1) is L1 + L2 + L3 - L4 = 0.
# This implies L4 = L1 + L2 + L3.
residue_eq = L1 + L2 + L3
print(f"The residue condition gives L4 = {residue_eq}.")

# Substitute this condition into the polynomial for Z_{0,4} to get Z_{0,3,1}.
Z_0_3_1 = Z_0_4.subs(L4, residue_eq)

# Expand the resulting polynomial to find its form.
Z_0_3_1_expanded = sympy.expand(Z_0_3_1)

print(f"\nSubstituting the residue condition, the polynomial Z_{{0,3,1}} becomes:")
# The prompt asks to output each number in the final equation.
# We will print the expanded polynomial, which shows all its terms.
# sympy's string representation uses ** for power.
final_poly_str = str(Z_0_3_1_expanded).replace('**', '^')
print(f"Z_{{0,3,1}} = {final_poly_str}")

# The degree of a multivariate polynomial is the maximum sum of exponents in any term.
# We can find this from the expanded form or use sympy's degree function.
final_degree = sympy.degree(Z_0_3_1_expanded, gen=[L1, L2, L3])
print(f"\nThe degree of the polynomial Z_{{0,3,1}} is {final_degree}.")
print("-" * 20)
