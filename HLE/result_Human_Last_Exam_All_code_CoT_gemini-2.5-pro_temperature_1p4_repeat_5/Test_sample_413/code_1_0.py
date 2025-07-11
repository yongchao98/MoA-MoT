from fractions import Fraction

# From our derivation, we found the function f(x) has the form (x-n)^3 - (x-n)
# with n = 1/4. We need to compute f(3).
n = Fraction(1, 4)
x = Fraction(3, 1)

# The expression for f(3) is (3 - 1/4)^3 - (3 - 1/4)
y = x - n

# Calculate the terms of the equation
# y = 3 - 1/4 = 11/4
term1 = y**3
term2 = y

# Calculate the final result
result = term1 - term2

# For printing the equation, we need a common denominator for the subtraction.
# The common denominator is 64.
term2_num_common = term2.numerator * (result.denominator // term2.denominator)
term2_den_common = result.denominator

# Output the equation and the final result as a fraction.
print(f"f(3) = (3 - 1/4)^3 - (3 - 1/4)")
print(f"f(3) = ({y.numerator}/{y.denominator})^3 - ({y.numerator}/{y.denominator})")
print(f"f(3) = {term1.numerator}/{term1.denominator} - {term2.numerator}/{term2.denominator}")
print(f"f(3) = {term1.numerator}/{term1.denominator} - {term2_num_common}/{term2_den_common}")
print(f"f(3) = {result.numerator}/{result.denominator}")