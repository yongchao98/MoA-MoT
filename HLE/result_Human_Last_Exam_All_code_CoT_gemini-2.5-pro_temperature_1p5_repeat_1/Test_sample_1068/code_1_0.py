import math

# The problem asks for the value of l(14).
# Our derived formula is l(p) = ln(Gamma(p+0.5)) - ln(Gamma(p)) + ln(8/sqrt(pi)).
# We need to calculate this for p = 14.
p = 14

# The ln(Gamma(x)) function is available as math.lgamma(x).
# So, l(14) = lgamma(14.5) - lgamma(14) + log(8/sqrt(pi))
# which is lgamma(14.5) - lgamma(14) + log(8) - 0.5*log(pi)

# Let's verify the derivation by transforming the expression algebraically
# to see if it simplifies to something that can be computed without lgamma.
# C = ln(8/sqrt(pi))
# l(14) = ln(Gamma(14.5)) - ln(Gamma(14)) + C
# Using Gamma(z+1/2) = (2z)! / (4^z * z!) * sqrt(pi)
# Let z=14. Gamma(14.5) = 28! / (4^14 * 14!) * sqrt(pi)
# l(14) = ln[ 28!*sqrt(pi)/(4^14*14!) ] - ln(13!) + ln(8/sqrt(pi))
# l(14) = ln(28!) + 0.5*ln(pi) - 14*ln(4) - ln(14!) - ln(13!) + ln(8) - 0.5*ln(pi)
# l(14) = ln(28!) - 28*ln(2) - ln(14*13!) - ln(13!) + 3*ln(2)
# l(14) = ln(28!) - ln(14) - 2*ln(13!) - 25*ln(2)
# This is ln(Gamma(29)) - ln(14) - 2*ln(Gamma(14)) - 25*ln(2)
# Another path: l(p) = ln(Gamma(p+1/2)) - ln(Gamma(p)) + ln(8/sqrt(pi)) can be simplified to:
# l(14) = ln( (27!) / (2**24 * (13!)**2) )
# This is lgamma(28) - 24*log(2) - 2*lgamma(14)

term1 = math.lgamma(28)
term2 = 24 * math.log(2)
term3 = 2 * math.lgamma(14)

result = term1 - term2 - term3

# The final result is the numerical value.
# The question also asks to show the final equation, which involves the numbers used.
p_val = 14
# The final expression for l(14) is ln(Gamma(27+1)) - 24*ln(2) - 2*ln(Gamma(13+1))
# We print the calculation with the final values.
print(f"The calculation is math.lgamma({p_val*2}) - {24}*math.log(2) - {2}*math.lgamma({p_val})")
print(f"Value of math.lgamma(28): {term1}")
print(f"Value of 24*math.log(2): {term2}")
print(f"Value of 2*math.lgamma(14): {term3}")
print(f"Final value of l(14) = {term1} - {term2} - {term3}")
print(f"l(14) = {result}")
print("<<<" + str(result) + ">>>")