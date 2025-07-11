import math

# The problem simplifies to the evaluation of the expression:
# (3/2) * 10^(10/3) + 37/4
# We will calculate this and show each component as requested.

# Define the numbers in the final equation
term1_coeff_num = 3
term1_coeff_den = 2
term1_base = 10
term1_exp_num = 10
term1_exp_den = 3

term2_num = 37
term2_den = 4

# Calculate the values of the components
term1_coeff = term1_coeff_num / term1_coeff_den
term1_exp = term1_exp_num / term1_exp_den
term1_val = term1_coeff * (term1_base ** term1_exp)

term2_val = term2_num / term2_den

# Calculate the final result
final_result = term1_val + term2_val

# Print the final equation with all its numeric components
print("The final expression is: (A/B) * C^(D/E) + F/G")
print(f"The numbers in this equation are: A={term1_coeff_num}, B={term1_coeff_den}, C={term1_base}, D={term1_exp_num}, E={term1_exp_den}, F={term2_num}, G={term2_den}\n")

print("The final calculation is:")
print(f"({term1_coeff_num}/{term1_coeff_den}) * {term1_base}^({term1_exp_num}/{term1_exp_den}) + ({term2_num}/{term2_den})")
print(f"= {term1_val:.4f} + {term2_val:.4f}")
print(f"= {final_result:.4f}")