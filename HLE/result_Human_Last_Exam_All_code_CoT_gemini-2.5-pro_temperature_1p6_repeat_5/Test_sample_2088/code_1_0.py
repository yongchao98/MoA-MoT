# The value of the definite integral was determined analytically.
# This script presents the final calculation based on that result.

term_base = 12
term_power = 4
integral_result_num_base = 5
integral_result_den = 12
final_answer = 5

# We need to print the final equation showing how the terms cancel out.
# The equation is: (12)^4 * ( (5^(1/4)) / 12 )^4 = 5
print(f"The final calculation is:")
print(f"({term_base})^{term_power} * (({integral_result_num_base}^(1/4) / {integral_result_den}))^{term_power} = {final_answer}")