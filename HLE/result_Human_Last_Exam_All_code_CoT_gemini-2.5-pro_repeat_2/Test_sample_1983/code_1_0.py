import math

# Step 1: Define the constants from the problem statement.
# We interpret "lambda_2 = 0.5 * lambda_1 = 0.5" as lambda_1 = 0.5 and lambda_2 = 0.5,
# as other interpretations lead to the inapplicability of the theorems in the reference.
k1_exp = 3000  # k1 = 10^3000
k2_exp = 500   # k2 = 10^500
lambda1 = 0.5
lambda2 = 0.5
h_norm = 1000.0

# Step 2: Calculate the limit superior and limit inferior of the norm of the solution x_n.
# Based on Theorem 4.2 from the reference [1], we use the sharp bounds.
# L_plus_inf = lim_sup_{n->+inf} ||x_n|| = (k1 / (1 - lambda1)) * |||h|||
# We will use logarithms to handle the large numbers k1 and k2.
log10_L_plus_inf = k1_exp + math.log10(h_norm / (1 - lambda1))

# L_minus_inf = lim_inf_{n->-inf} ||x_n|| = (k2 / (1 - lambda2)) * |||h|||
log10_L_minus_inf = k2_exp + math.log10(h_norm / (1 - lambda2))

# Step 3: Calculate the two terms of the final expression.
# The expression is E = 100 * log10( (1/3) * L_plus_inf ) + 10 * log10( (1/3) * L_minus_inf )

# First term: 100 * log10( (1/3) * L_plus_inf )
# This is 100 * (log10(1/3) + log10(L_plus_inf))
term1 = 100 * (math.log10(1/3) + log10_L_plus_inf)

# Second term: 10 * log10( (1/3) * L_minus_inf )
# This is 10 * (log10(1/3) + log10(L_minus_inf))
term2 = 10 * (math.log10(1/3) + log10_L_minus_inf)

# Step 4: Calculate the final result.
final_result = term1 + term2

# Step 5: Print the breakdown of the calculation as requested.
print("Problem to solve: 100 * limsup log10(1/3 ||x_n||) + 10 * liminf log10(1/3 ||x_n||)")
print("-" * 20)
print(f"Calculated value for the first term: 100 * log10(1/3 * limsup ||x_n||) = {term1}")
print(f"Calculated value for the second term: 10 * log10(1/3 * liminf ||x_n||) = {term2}")
print("-" * 20)
print(f"The final equation is: {term1} + {term2} = {final_result}")
print("-" * 20)
print(f"Final Result: {final_result}")

# Let's show the full analytic calculation as well
# Term 1 = 100 * (log10(1/3) + log10(10^3000 * 1000 / 0.5)) = 100 * (log10(1/3) + log10(2 * 10^3003)) = 100 * (log10(2/3) + 3003)
# Term 2 = 10 * (log10(1/3) + log10(10^500 * 1000 / 0.5)) = 10 * (log10(1/3) + log10(2 * 10^503)) = 10 * (log10(2/3) + 503)
# Total = 100*log10(2/3) + 300300 + 10*log10(2/3) + 5030 = 110*log10(2/3) + 305330
analytic_result = 110 * math.log10(2.0/3.0) + 305330
# This is just to verify the code's logic. The output above is self-contained.
# print(f"Analytic check: {analytic_result}")