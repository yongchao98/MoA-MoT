import scipy.stats
import numpy as np

# --- Problem Parameters ---
# We have 1 instrument.
num_instruments = 1  # K

# We want the relative asymptotic bias to be less than 10%.
relative_bias_max = 0.10  # r

# We want 95% confidence.
confidence_level = 0.95
quantile = confidence_level

# --- Step 1: Explain the governing equation ---
print("The relationship between maximum relative TSLS bias (r), the number of instruments (K), and the non-centrality parameter (λ) of the first-stage F-statistic is:")
print("r = 1 / (1 + λ / K)")
print("-" * 50)

# --- Step 2: Calculate the required non-centrality parameter (λ) ---
# We rearrange the formula to solve for λ at the boundary condition.
# λ = K * (1/r - 1)
ncp = num_instruments * (1 / relative_bias_max - 1)

print(f"To ensure the relative bias is less than {relative_bias_max}, we solve for the boundary λ:")
print(f"λ = {num_instruments} * (1 / {relative_bias_max} - 1)")
print(f"λ = {num_instruments} * ({1/relative_bias_max} - 1)")
print(f"λ = {ncp:.2f}")
print("-" * 50)

# --- Step 3: Find the F-statistic critical value ---
# We need to find the critical value from the F-distribution that provides 95% confidence.
# This corresponds to the 95th percentile of the non-central F-distribution.
dfn = num_instruments
dfd = 1e9  # Using a large number to approximate infinity for the asymptotic case.

print(f"To have {confidence_level*100}% confidence that the instrument is strong enough (i.e., true λ > {ncp:.2f}),")
print(f"our observed F-statistic must exceed the {int(quantile*100)}th percentile of the distribution under the null hypothesis.")
print("\nThis distribution is the non-central F-distribution with parameters:")
print(f"  Numerator degrees of freedom (K) = {dfn}")
print(f"  Denominator degrees of freedom      = ∞ (approximated by {dfd})")
print(f"  Non-centrality parameter (λ)      = {ncp:.2f}")
print("-" * 50)

# --- Step 4: Perform the calculation and output the result ---
# Use the Percent Point Function (PPF), which is the inverse of the CDF.
critical_value = scipy.stats.ncf.ppf(quantile, dfn, dfd, ncp)
final_answer = round(critical_value)

print(f"The calculated critical value for the F-statistic is: {critical_value:.4f}")
print("\nThis is the well-known Stock-Yogo critical value for this case.")
print(f"Rounding to the nearest whole number, the minimum F-statistic required is: {int(final_answer)}")
