import math

# The coefficients for the correction terms in the logarithm's expansion.
# These come from the Euler-Maclaurin formula.
# The term for n^-2 comes from B_4 / (4!) * f'''(n)
# The term for n^-4 comes from B_6 / (6!) * f^(5)(n)

# First term coefficient: 1/720
c1_num = 1
c1_den = 720

# Second term coefficient: -1/5040
c2_num = 1
c2_den = 5040

# Construct and print the formula for P(n)
# The formula is constructed to meet the O(n^-6) relative error requirement.
# We use the exponential form as it is more compact and directly represents
# the correction derived from the logarithmic expansion.
formula = f"P(n) = exp( {c1_num} / ({c1_den}*n**2) - {c2_num} / ({c2_den}*n**4) )"

print(formula)