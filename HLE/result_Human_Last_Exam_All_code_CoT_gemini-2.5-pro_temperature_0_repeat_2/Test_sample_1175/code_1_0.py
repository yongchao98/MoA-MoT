import math

# After solving the differential equation and simplifying the expression,
# the final calculation is (3/2) * 10^(10/3) + 37/4.

# Define the components of the expression using their fractional form for clarity.
coeff1_num, coeff1_den = 3, 2
base1 = 10
exp1_num, exp1_den = 10, 3
term2_num, term2_den = 37, 4

# Calculate the floating-point values for computation.
coeff1 = coeff1_num / coeff1_den
exp1 = exp1_num / exp1_den
term2 = term2_num / term2_den

# Calculate the value of the first term.
term1_val = coeff1 * (base1 ** exp1)

# Calculate the final result by adding the two terms.
result = term1_val + term2

# Print the final equation with the numbers and the result.
print(f"The final expression to evaluate is: ({coeff1_num}/{coeff1_den}) * {base1}^({exp1_num}/{exp1_den}) + ({term2_num}/{term2_den})")
print(f"Value of the first term: {term1_val}")
print(f"Value of the second term: {term2}")
print(f"Final result: {result}")