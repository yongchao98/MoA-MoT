import math

# Calculate the value of the limit expression
# alpha = (13 + sqrt(165))/2, ln(alpha) = acosh(13/2) = acosh(6.5)
# gamma = (5 + sqrt(21))/2, ln(gamma) = acosh(5/2) = acosh(2.5)
limit_coeff = 2 / math.acosh(6.5) + 2 / math.acosh(2.5)

# The problem asks for the integer part of 10^4 * limit_coeff
result = 10**4 * limit_coeff

print(f"The value of alpha is (13+sqrt(165))/2")
print(f"The value of gamma is (5+sqrt(21))/2")
print(f"The limit L = F(N)/ln(N) is 2/ln(alpha) + 2/ln(gamma)")
print(f"This is equal to 2/arccosh(6.5) + 2/arccosh(2.5)")
print(f"Numerically, L is approximately {limit_coeff}")
print(f"The value 10^4 * L is approximately {result}")
print(f"The integer part of 10^4 * L is {math.floor(result)}")