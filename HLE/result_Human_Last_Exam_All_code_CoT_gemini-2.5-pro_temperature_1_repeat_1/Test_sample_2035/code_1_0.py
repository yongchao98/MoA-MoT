import math

# Step 1: Define the parameters from the problem description.
n = 4048

# Step 2: Calculate the sum of alpha_i^2.
# We are given alpha_i = sqrt(i), so alpha_i^2 = i.
# We need to calculate the sum of i for i from 1 to n.
# This is the sum of the first n integers, which has a well-known formula: n * (n + 1) / 2.
# We use integer division // to ensure the result is an integer.
sum_alpha_sq = n * (n + 1) // 2

# Step 3: Calculate the first term of the final expression.
# Based on the derivation, the expression simplifies to (n/2) * sum(alpha_i^2).
# Let's call this term A.
term_A = (n // 2) * sum_alpha_sq

# Step 4: Define the second term of the final expression.
# This is given as 10^15. Let's call this term B.
term_B = 10**15

# Step 5: Calculate the final result by summing the two terms.
final_result = term_A + term_B

# Step 6: Output the components of the final equation and the result, as requested.
# The simplified equation is A + B = final_result.
print(f"The simplified expression to be calculated is of the form A + B.")
print(f"The value of term A, derived from (n/2) * sum(alpha_i^2), is: {term_A}")
print(f"The value of term B is: {term_B}")
print(f"The final equation is: {term_A} + {term_B} = {final_result}")
