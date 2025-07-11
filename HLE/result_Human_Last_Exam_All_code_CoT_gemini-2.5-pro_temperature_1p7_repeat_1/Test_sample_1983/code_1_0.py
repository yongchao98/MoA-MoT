import math

# Step 1 & 3: Define parameters based on our interpretation
k1_exp = 3000
k2_exp = 500
h_norm = 1000

# We assume lambda_1 = 0.5 based on the analysis that lambda_1=1 leads to a divergent result.
# This implies lambda_2 = 0.5 * lambda_1.
lambda1 = 0.5
lambda2 = 0.5 * lambda1

# Step 2 & 4: Calculate the arguments for the logarithms.
# The limits are too large to be represented as standard floats. We compute their logarithms directly.
# Let L_pos_inf = (k1 / (1 - lambda1)) * h_norm
# log10(L_pos_inf / 3) = log10( (10^k1_exp / (1 - lambda1) * h_norm) / 3)
# = log10(10^3000 / 0.5 * 1000 / 3) = log10(2 * 10^3003 / 3) = log10(2/3) + 3003
log_val1 = math.log10(2.0 / 3.0) + 3003

# Let L_neg_inf = (k2 / (1 - lambda2)) * h_norm
# log10(L_neg_inf / 3) = log10( (10^k2_exp / (1 - lambda2) * h_norm) / 3)
# = log10(10^500 / 0.75 * 1000 / 3) = log10( (4/3) * 10^503 / 3) = log10(4/9) + 503
log_val2 = math.log10(4.0 / 9.0) + 503

# Step 5: Compute the final expression as per the problem statement
# Expression E = 100 * log_val1 + 10 * log_val2
result = 100 * log_val1 + 10 * log_val2

# Print the final equation with the computed numerical values
print(f"100 * ({log_val1}) + 10 * ({log_val2}) = {result}")
<<<305308.8690489126>>>