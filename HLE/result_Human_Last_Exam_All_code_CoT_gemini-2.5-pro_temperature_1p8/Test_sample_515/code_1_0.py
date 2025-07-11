import math

# Define the constants
L = 3600.0
# Euler-Mascheroni constant
gamma = 0.5772156649
# Constant C in the asymptotic expansion of the potential kernel
C = (2 * gamma + math.log(8)) / math.pi
# pi
pi = math.pi

# Numerator of the probability formula is approximately C
numerator = C

# Denominator of the probability formula
denominator_term1 = (4 / pi) * math.log(L)
denominator_term2 = 2 * C
denominator = denominator_term1 + denominator_term2

# The probability P
P = numerator / denominator

print(f"The starting point is (0,1). The target z_A is (3600,0).")
print(f"The probability is given by the formula: C / ( (4/pi)*ln(L) + 2*C )")
print(f"L = {L}")
print(f"pi = {pi}")
print(f"C = (2*gamma + ln(8))/pi = {C}")
print(f"Numerator = C = {numerator}")
print(f"Denominator = (4/pi)*ln({L}) + 2*C = {denominator_term1} + {denominator_term2} = {denominator}")
print(f"Probability = {numerator} / {denominator} = {P}")
print(f"The approximate answer with two significant digits is: {P:.2g}")