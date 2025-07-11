import sympy

# Define m as a symbolic variable to represent it in the formula
m = sympy.Symbol('m')

# Number of successful pairs (i, j) based on our combinatorial model.
# S_m = C(m+2, 2)
num_successful_pairs = (m + 1) * (m + 2) / 2

# Total number of pairs (i, j) that can be chosen.
# T_m = C(4m+2, 2)
total_pairs = (4*m + 2) * (4*m + 1) / 2

# The probability P_m is the ratio S_m / T_m.
# We can simplify the expression.
P_m_numerator = (m + 1) * (m + 2)
P_m_denominator = 2 * (2*m + 1) * (4*m + 1)

# Display the final formula for P_m
# As requested, we show the numbers and structure of the final formula.
print("The probability P_m for a given positive integer m is calculated as:")
print(f"P_m = ( (m + 1) * (m + 2) ) / ( 2 * (2*m + 1) * (4*m + 1) )")