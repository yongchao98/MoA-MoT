import sympy

# Define symbols
n = sympy.Symbol('n')
L = sympy.Symbol('L')

# Coefficients for the terms in P(n)
# For the n^-2 term: (c2_L2*L**2 + c2_L1*L + c2_L0) / c2_den
c2_L2 = 3
c2_L1 = -2
c2_L0 = 2
c2_den = 24

# For the n^-3 term: (c3_L3*L**3 + c3_L2*L**2 + c3_L1*L) / c3_den
c3_L3 = 1
c3_L2 = -2
c3_L1 = 6
c3_den = 48

# Construct the expression for P(n)
term1 = (c2_L2*L**2 + c2_L1*L + c2_L0) / (c2_den * n**2)
term2 = (c3_L3*L**3 + c3_L2*L**2 + c3_L1*L) / (c3_den * n**3)
P_n = term1 + term2

# Print the formula in a readable format
print("The formula for P(n) is:")
# We manually construct the string to ensure all numbers are explicitly shown
# as requested.
formula_str = (
    f"({c2_L2}*L^2 - {abs(c2_L1)}*L + {c2_L0})/({c2_den}*n^2) + "
    f"({c3_L3}*L^3 - {abs(c3_L2)}*L^2 + {c3_L1}*L)/({c3_den}*n^3)"
)
# Adjust signs for pretty printing
formula_str = formula_str.replace("- -", "+ ").replace("+ -", "- ")
print(f"P(n) = {formula_str}")

# Final answer in specified format
final_answer_formula = "(3*L^2 - 2*L + 2)/(24*n^2) + (L^3 - 2*L^2 + 6*L)/(48*n^3)"
# The prompt is a bit ambiguous whether the output should be just the formula
# or a full python block. The block above shows the derivation and prints the formula.
# The final extracted answer is just the formula itself.
