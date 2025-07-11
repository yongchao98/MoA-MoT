import scipy.stats
import math

#
# This script calculates the minimum F-statistic required to be 95% confident
# that the relative bias of a TSLS estimate is less than 10%, given one
# instrumental variable and one endogenous variable.
# This value is known as the Stock-Yogo critical value for this scenario.
#

# --- Step 1: Define the parameters for the calculation ---

# The confidence level is 95%, so we will find the 0.95 quantile of the distribution.
q = 0.95

# The degrees of freedom (df) for the asymptotic non-central chi-squared distribution
# is equal to the number of instrumental variables.
df = 1

# The requirement is that the relative bias is less than 10% (0.10).
# The relative bias is approximated by 1 / E[F].
# So, we need E[F] > 1 / 0.10 = 10.
# The null hypothesis for our test is E[F] <= 10.
#
# The expected F-statistic E[F] is related to the non-centrality parameter (nc)
# of its asymptotic distribution by E[F] ≈ 1 + nc/df.
# From this, the non-centrality parameter under the null hypothesis (E[F]=10) is:
# nc = (10 - 1) * df = 9.
nc = 9

# --- Step 2: Calculate the critical F-statistic ---

# The minimum required F-statistic is the critical value from the non-central
# chi-squared distribution with the parameters defined above. We use the
# Percent Point Function (PPF), which is the inverse of the CDF, to find this value.
f_statistic_critical = scipy.stats.ncx2.ppf(q, df, nc)

# The question asks for the nearest whole number.
f_statistic_rounded = round(f_statistic_critical)

# --- Step 3: Print the final equation and the result ---

print("This is the derivation for the minimum required F-statistic.")
print("The 'final equation' shows the statistical calculation and its parameters:")
print(f"Minimum F-statistic = round( P.P.F. of non-central chi-squared(q={q}, df={df}, nc={nc}) )")
print("\nPlugging in the numbers and the result from the calculation:")
print(f"{f_statistic_rounded} = round( {f_statistic_critical:.4f} )")