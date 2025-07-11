# Based on the analysis of the graph's features, we construct the equation for f(x).
# A logical inconsistency in the provided graph is resolved by assuming a minor plotting error.
# The final equation represents all key features: roots, vertical asymptotes, and the slant asymptote.

# The equation is formed by terms with specific exponents (powers).
# The prompt requests to output each number (exponent) in the final equation.

# Numerator Term 1 (from root at x=d): (x - d) with an exponent of 2
power_d = 2

# Numerator Term 2 (from roots at x=b, x=-b): (x^2 - b^2) with an exponent of 1
# This itself contains terms x^2 and b^2
power_b_group = 1
power_x_in_b = 2
power_b_in_b = 2

# Denominator Term 1 (from asymptote at x=a): (x - a) with an exponent of 1
power_a = 1

# Denominator Term 2 (from asymptote at x=c): (x - c) with an exponent of 2
power_c = 2

# Print the final equation using f-strings to insert the exponents.
print("f(x) = ( (x - d)^{} * (x^{} - b^{})^{} ) / ( (x - a)^{} * (x - c)^{} )".format(
    power_d, power_x_in_b, power_b_in_b, power_b_group, power_a, power_c
))