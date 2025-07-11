import math

# Final expression to calculate is (3/2) * 10^(10/3) + 37/4

# Define the terms of the equation
term1_coeff_num = 3
term1_coeff_den = 2
term1_base = 10
term1_exp_num = 10
term1_exp_den = 3

term2_num = 37
term2_den = 4

# Calculate the value of each term
term1_value = (term1_coeff_num / term1_coeff_den) * (term1_base**(term1_exp_num / term1_exp_den))
term2_value = term2_num / term2_den

# Calculate the final result
result = term1_value + term2_value

# Print the final equation with the numerical value of each term
# as requested.
print("The simplified expression is: (3/2) * 10^(10/3) + 37/4")
print(f"The first term is ({term1_coeff_num}/{term1_coeff_den}) * {term1_base}^({term1_exp_num}/{term1_exp_den}) = {term1_value}")
print(f"The second term is {term2_num}/{term2_den} = {term2_value}")
print("\nThe final equation is:")
print(f"{term1_value} + {term2_value} = {result}")
