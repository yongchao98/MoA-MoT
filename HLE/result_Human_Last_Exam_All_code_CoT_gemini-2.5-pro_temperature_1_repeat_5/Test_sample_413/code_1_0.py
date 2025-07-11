from fractions import Fraction

# From the analysis, the roots of the polynomial f(x) are of the form
# {n, n, s} where n is an integer.
# The condition f'(-1/4) = -1/4 gives a relation between the coefficients a and b,
# which in turn relates n and s.
# After testing integer values for n, we found a unique solution that satisfies all
# the problem's conditions: n=1 and s = -39/40.
# The roots are therefore {1, 1, -39/40}.
# The polynomial is f(x) = (x-1)^2 * (x - (-39/40)).
# Let's find the coefficients a, b, c to form the polynomial f(x) = x^3 + ax^2 + bx + c.
# a = -(s+2n) = -(-39/40 + 2) = -41/40
# b = n^2+2ns = 1 + 2*(-39/40) = 1 - 39/20 = -19/20
# c = -n^2*s = -1*(-39/40) = 39/40
a = Fraction(-41, 40)
b = Fraction(-19, 20)
c = Fraction(39, 40)

# We need to compute f(3).
x = 3
f_3 = x**3 + a * x**2 + b * x + c

# The final equation is f(3) = 3^3 + (-41/40)*3^2 + (-19/20)*3 + 39/40
# We want to print each number in the equation.
print(f"f(3) = {x}^3 + ({a.numerator}/{a.denominator})*{x}^2 + ({b.numerator}/{b.denominator})*{x} + {c.numerator}/{c.denominator}")

# Now we calculate and print the result.
# f(3) = 27 - (41*9)/40 - (19*3)/20 + 39/40
# f(3) = 27 - 369/40 - 57/20 + 39/40
# To sum them up, we find a common denominator, which is 40.
# f(3) = (27*40)/40 - 369/40 - (57*2)/40 + 39/40
# f(3) = (1080 - 369 - 114 + 39) / 40
# f(3) = (1119 - 483) / 40
# f(3) = 636 / 40
# f(3) = 159 / 10
print(f"f(3) = {f_3.numerator}/{f_3.denominator}")

final_answer = f_3
#<<<159/10>>>