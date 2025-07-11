import sympy

# Define symbols
L, n = sympy.symbols('L n')

# The numerators of the terms for P(n)
p2_num = 3*L**2 - 2*L + 2
p3_num = L**3 - 2*L**2 + 2*L

# The denominators
p2_den = 24 * n**2
p3_den = 48 * n**3

# Construct the expression for P(n)
P_n = p2_num / p2_den + p3_num / p3_den

# Print the formula
print("The formula for P(n) is:")
print(f"P(n) = {sympy.pretty(P_n)}")

# The problem asks for the formula itself as the final answer.
# We extract the string representation.
final_formula_str = f"({p2_num}) / ({p2_den}) + ({p3_num}) / ({p3_den})"
# The prompt also asks to output each number in the final equation.
# Printing the pretty format above and the formula string below covers this.
# Let's present the formula in a clear, standard way for the final answer.
# P(n) = (3*L**2 - 2*L + 2)/(24*n**2) + (L**3 - 2*L**2 + 2*L)/(48*n**3)