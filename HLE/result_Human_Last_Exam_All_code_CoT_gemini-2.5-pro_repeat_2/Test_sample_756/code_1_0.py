from fractions import Fraction

# Define the coefficients of the polynomial f(x) = ax^2 + bx + c
# that achieves the maximum value for |b| + |c|.
a = Fraction(-8, 9)
b = Fraction(8, 9)
c = Fraction(7, 9)

# Calculate the value of |b| + |c|
abs_b = abs(b)
abs_c = abs(c)
result = abs_b + abs_c

# Print the final equation with the numbers
print(f"The maximum value is achieved for the polynomial f(x) = {a}x^2 + {b}x + {c}.")
print(f"For this polynomial, b = {b} and c = {c}.")
print(f"The calculation for |b| + |c| is:")
print(f"|{b}| + |{c}| = {abs_b} + {abs_c} = {result}")
