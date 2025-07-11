from fractions import Fraction

# From the analysis, the coefficients are determined to be:
a = Fraction(-3, 2)
b = Fraction(-19, 16)
c = Fraction(21, 16)

# The function is f(x) = x^3 + ax^2 + bx + c
# We need to compute f(3)
x = 3
f_3 = x**3 + a * x**2 + b * x + c

# The final equation is f(3) = 3^3 + (-3/2)*3^2 + (-19/16)*3 + 21/16
# Print out each part of the calculation
p1 = x**3
p2 = a * x**2
p3 = b * x
p4 = c

# To avoid float issues and present the exact fraction
p1_frac = Fraction(p1)
p2_frac = a * Fraction(x**2)
p3_frac = b * Fraction(x)
p4_frac = c

print(f"The coefficients are a = {a}, b = {b}, c = {c}.")
print(f"The function is f(x) = x^3 + ({a})x^2 + ({b})x + ({c}).")
print("We want to compute f(3).")
print(f"f(3) = 3^3 + ({a})*3^2 + ({b})*3 + ({c})")
print(f"f(3) = {p1_frac} + ({p2_frac}) + ({p3_frac}) + ({p4_frac})")
print(f"f(3) = {p1_frac + p2_frac + p3_frac + p4_frac}")
print(f"The exact value of f(3) is {f_3}.")

final_answer = f_3
# The final result is a fraction.
# Let's output it in the specified format.
# No, I should not do this. I need to output the answer as a number or fraction string.
# <<<45/4>>>