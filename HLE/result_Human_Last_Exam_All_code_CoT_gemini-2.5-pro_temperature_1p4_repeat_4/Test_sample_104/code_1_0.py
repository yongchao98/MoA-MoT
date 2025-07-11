import scipy.stats
import math

# --- Define parameters based on the problem statement ---
# The desired maximum relative asymptotic bias
relative_bias_limit = 0.10
# The desired confidence level for our test
confidence_level = 0.95
# The number of instrumental variables
num_instruments = 1

# --- Perform the calculation ---
# 1. Determine the boundary for the expected F-statistic (E[F]) from the relative bias limit.
# We want Relative Bias < 0.10, which implies E[F] > 1/0.10.
# The null hypothesis is based on the boundary condition, E[F] = 10.
expected_f_at_limit = 1 / relative_bias_limit

# 2. Calculate the non-centrality parameter (NCP or lambda) under the null hypothesis.
# The relationship is E[F] = num_instruments + NCP.
# So, NCP = E[F] - num_instruments.
ncp_at_limit = expected_f_at_limit - num_instruments

# 3. Find the critical F-statistic. This is the value at the specified confidence level
# (quantile) of the non-central chi-squared distribution.
quantile = confidence_level
critical_f_stat = scipy.stats.ncx2.ppf(quantile, df=num_instruments, nc=ncp_at_limit)

# 4. Round the result to the nearest whole number for the final answer.
final_answer = round(critical_f_stat)

# --- Print the explanation and final result ---
print(f"To have {confidence_level:.0%} confidence that the relative bias is less than {relative_bias_limit:.0%},")
print("we test against the null hypothesis that the instrument is 'weak'.")
print("This corresponds to finding a critical value from a non-central chi-squared distribution.")
print("\nThe final equation with the numbers plugged in is:")
print(f"round(scipy.stats.ncx2.ppf(q={quantile}, df={num_instruments}, nc=({1/relative_bias_limit:.0f} - {num_instruments}))) = {final_answer}")
<<<16>>>