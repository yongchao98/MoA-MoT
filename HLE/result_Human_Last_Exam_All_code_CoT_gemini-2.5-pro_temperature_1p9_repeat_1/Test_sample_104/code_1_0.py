import scipy.stats

# Define the parameters for the problem
confidence_level = 0.95
max_relative_bias = 0.10
number_of_instruments = 1

# Step 1: From the relative bias condition, determine the boundary for the
# expected F-statistic. We want Bias < 0.10, which implies E[F] > 1/0.10.
boundary_expected_F = 1 / max_relative_bias

# Step 2: Determine the parameters of the asymptotic distribution (non-central chi-squared)
# under the boundary condition E[F] = 10.
# Degrees of freedom (df) = number of instruments
df = number_of_instruments
# Non-centrality parameter (ncp) = E[F] - df
ncp = boundary_expected_F - df

# Step 3: Calculate the critical value. This is the 95th percentile of the
# chi-squared distribution with df=1 and ncp=9. If we observe an F-statistic
# greater than this value, we can be 95% confident that the true E[F] > 10.
F_statistic_critical = scipy.stats.ncx2.ppf(confidence_level, df, ncp)

# Step 4: Round the result to the nearest whole number as requested.
final_answer = round(F_statistic_critical)

# Step 5: Output the numbers used in the final calculation and the answer.
print("To find the minimum F-statistic, we find the 95th percentile of a non-central chi-squared distribution.")
print("The final calculation is F_statistic = ppf(q, df, nc) where:")
print(f"q (confidence level) = {confidence_level}")
print(f"df (degrees of freedom) = {df}")
print(f"nc (non-centrality parameter) = {ncp}")
print(f"\nThe calculated critical F-statistic is {F_statistic_critical:.3f}.")
print(f"The minimum F-statistic required, rounded to the nearest whole number, is: {final_answer}")