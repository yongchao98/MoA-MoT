import math

# Step 1: Define the given parameters.
n = 4048
power_of_ten = 15

# Step 2: Calculate the sum of alpha_i^2.
# As derived in the plan, the problem simplifies to calculating the sum of alpha_i^2 + 10^15.
# Given alpha_i = sqrt(i), alpha_i^2 = i.
# We need to calculate the sum of integers from 1 to n.
# The formula for the sum of the first n integers is n * (n + 1) / 2.
sum_alpha_sq = n * (n + 1) // 2

# Step 3: Calculate the constant term.
constant_term = 10**power_of_ten

# Step 4: Calculate the final result.
final_result = sum_alpha_sq + constant_term

# Step 5: Print the components of the final equation and the result.
print(f"The simplified final equation is: (sum of alpha_i^2 from i=1 to n) + 10^{power_of_ten}")
print(f"The value for 'sum of alpha_i^2' (where n={n}) is: {sum_alpha_sq}")
print(f"The value for '10^{power_of_ten}' is: {constant_term}")
print(f"The final result is: {final_result}")
