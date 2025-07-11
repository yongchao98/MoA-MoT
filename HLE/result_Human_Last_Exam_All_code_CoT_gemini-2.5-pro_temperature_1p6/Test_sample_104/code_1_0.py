import scipy.stats

# --- Plan Summary ---
# We want to find the minimum F-statistic that provides 95% confidence that the
# relative asymptotic bias of the TSLS estimator is less than 10%.
#
# 1. The relative bias is approximated by 1 / F_pop, where F_pop is the
#    population F-statistic.
# 2. For a bias < 10% (0.1), we require F_pop > 1 / 0.1 = 10.
# 3. F_pop is equal to the non-centrality parameter (NCP or λ) for a single
#    instrument (K=1). So, we require λ > 10.
# 4. We set up a hypothesis test where the null hypothesis is H0: λ <= 10.
# 5. The critical value for the F-statistic is the 95th percentile of the
#    F-statistic's distribution when λ = 10. Asymptotically, this distribution
#    is a non-central chi-squared distribution with df=K=1 and nc=λ=10.

# --- Calculation ---

# Define the parameters for the calculation based on the problem statement
confidence_level = 0.95
relative_bias_threshold = 0.10
num_instruments = 1

# From the approximation, calculate the required non-centrality parameter (NCP)
# under the null hypothesis boundary.
# ncp = num_instruments / relative_bias_threshold
ncp = num_instruments / relative_bias_threshold

# The degrees of freedom for the chi-squared distribution is the number of instruments
df = num_instruments

# Calculate the critical value (the F-statistic) using the Percent Point Function (PPF),
# which is the inverse of the CDF. We find the value for which the CDF is 0.95.
f_statistic_critical_value = scipy.stats.ncx2.ppf(confidence_level, df, ncp)

# Round the result to the nearest whole number as requested.
final_answer = round(f_statistic_critical_value)


print("--- Derivation of the Minimum F-Statistic ---")
print(f"Confidence Level = {confidence_level*100:.0f}%")
print(f"Maximum Relative Bias Threshold = {relative_bias_threshold*100:.0f}%")
print(f"Number of Instruments = {num_instruments}")
print("\n--- The Final 'Equation' ---")
print(f"To satisfy the condition, we need to be able to reject the null hypothesis that the Non-Centrality Parameter (λ) is less than or equal to {ncp:.0f}.")
print(f"The critical F-statistic is the {confidence_level*100:.0f}th percentile of the non-central chi-squared distribution with {df} degree(s) of freedom and a non-centrality parameter of {ncp:.0f}.")
print(f"Calculation: scipy.stats.ncx2.ppf({confidence_level}, df={df}, nc={ncp:.0f}) = {f_statistic_critical_value:.4f}")
print("\n--- Final Answer ---")
print("The minimum F-statistic required, rounded to the nearest whole number, is:")
print(int(final_answer))
