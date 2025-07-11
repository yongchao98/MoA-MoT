from fractions import Fraction

# Step 4: Define x0 and compute f(3)
x0 = Fraction(1, 4)
x = 3

# The function is f(x) = (x - x0)^3 - (x - x0)
val = Fraction(x) - x0
result = val**3 - val

# The problem asks for the final equation as well.
# f(3) = (3 - 1/4)^3 - (3 - 1/4)
# f(3) = (11/4)^3 - (11/4)
# f(3) = 1331/64 - 11/4
# f(3) = 1331/64 - 176/64
# f(3) = 1155/64

print("The function is f(x) = (x - 1/4)^3 - (x - 1/4).")
print("We want to compute f(3).")
print("f(3) = (3 - 1/4)^3 - (3 - 1/4)")
term1 = (Fraction(3) - Fraction(1, 4))
term1_cubed = term1**3
print(f"f(3) = ({term1})^3 - ({term1})")
print(f"f(3) = {term1_cubed} - {term1}")
term1_common_denom = term1.limit_denominator(64)
term1_cubed_common_denom = term1_cubed.limit_denominator(64)
print(f"f(3) = {term1_cubed_common_denom} - {term1_common_denom}")
final_numerator = term1_cubed.numerator * term1.denominator - term1.numerator * term1_cubed.denominator
final_denominator = term1_cubed.denominator * term1.denominator
# This is not quite right, need common denominator
final_numerator = term1_cubed.numerator - term1.numerator * 16
final_denominator = 64
print(f"f(3) = ({term1_cubed.numerator} - {term1.numerator*16}) / {final_denominator}")
print(f"f(3) = {final_numerator}/{final_denominator}")
print(f"The exact value of f(3) is {result}.")
